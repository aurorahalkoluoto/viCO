library(Seurat)
library(tidyverse)

obj <- readRDS("/scratch/project_2012655/pau/Organoids_MergedRawCounts.RDS")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

## Find highly variable features
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)


# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "/scratch/project_2012655/pau/nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Decision on whether to use ccdifference or no cell cycle regression at all in “ScaleData"
obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, features = VariableFeatures(object = obj))


# Examine and visualize PCA results a few different ways
print(obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(obj, dims = 1:2, reduction = "pca")

## PCA plot with cells colored by 10x library
dirProcessing <- "/scratch/project_2012655/Aurora/QC_plot01/"

plotToPdf <- DimPlot(obj, reduction = "pca", pt.size = .1, group.by = 'orig.ident')+ggtitle("Before Harmony")
pdf(file=paste0(dirProcessing,"pooled_beforeHarmony.pdf"))
plot(plotToPdf)
dev.off()

library(harmony)
library(tidyverse)

## define batch

obj$batch <- ifelse(grepl("Alz43cl2_", obj$orig.ident) | grepl("PLCG1_", obj$orig.ident), "B1", "B2")
obj$condition <- sapply(strsplit(obj$orig.ident, "_"), function(x) paste(x[-1], collapse="_"))
obj$donor <- sapply(strsplit(obj$orig.ident, "_"), function(x) paste(x[1], collapse="_"))


# obj <- obj %>% 
#   RunHarmony(c("orig.ident"), plot_convergence = TRUE)
set.seed(123)
obj <- obj %>% 
  RunHarmony(c("donor"), plot_convergence = TRUE)

harmony_embeddings <- Embeddings(obj, 'harmony')
harmony_embeddings[1:5, 1:5]

## PCA-corrected plot with cells colored by 10x library
plotToPdf <- DimPlot(obj, reduction = "pca", pt.size = .1, group.by = "donor") + ggtitle("Before Harmony")
plotToPdf2 <- DimPlot(obj, reduction = "harmony", pt.size = .1, group.by = "donor") + ggtitle("After Harmony")
plotToPdf3 <- DimPlot(obj, reduction = "harmony", pt.size = .1, group.by = "orig.ident") + ggtitle("After Harmony")+
  facet_wrap(~orig.ident, nrow=2)



pdf(file=paste0(dirProcessing,"pooled_AfterHarmony.pdf"))
plotToPdf+plotToPdf2
dev.off()

# Use pca-corrected embeddings from harmony to produce UMAP and run clustering
obj <- obj %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>%
  FindClusters(resolution = 0.1) %>%
  identity()
#was 0.5 and 0.1

## UMAP plot: Cells colored by unannotated  clusters
plotToPdf <- DimPlot(obj, reduction = "umap", label = FALSE, pt.size = .1, group.by="condition")+theme_gray()+
  theme(plot.title=element_text(hjust=0.5, size=13, face="bold"))+
  ggtitle("Condition")+theme_bw()+facet_wrap(~condition)

pdf(file=paste0(dirProcessing,"cluteringUMAP_notAnnotated_condition.pdf"))
plot(plotToPdf)
dev.off()

plotToPdf <- DimPlot(obj, reduction = "umap", label = FALSE, pt.size = .1, group.by="orig.ident")+theme_gray()+
  theme(plot.title=element_text(hjust=0.5, size=13, face="bold"))+
  ggtitle("UMAP per sample")+theme_bw()+facet_wrap(~orig.ident, ncol=3)

pdf(file=paste0(dirProcessing,"cluteringUMAP_notAnnotated_batch.pdf"))
plot(plotToPdf)
dev.off()


plotToPdf <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "batch",
  split.by = "condition",
  label = FALSE,
  pt.size = 0.1,
  ncol = 3
) +
  ggtitle("UMAP per sample") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        strip.text = element_text(size=12))

pdf(file=paste0(dirProcessing,"cluteringUMAP_notAnnotated_batch_condition.pdf"))
plot(plotToPdf)
dev.off()



#saveRDS(obj, file="/scratch/project_2012655/pau/Organoids_notAnnotated_donor.RDS")
saveRDS(obj, file="/scratch/project_2012655/Aurora/Organoids_notAnnotated_donor_01.RDS")


#########
#########


# MILOR
# R packages
.libPaths(c("/projappl/project_2012655/rpackages", .libPaths()))
libpath <- .libPaths()[1]

library(SingleCellExperiment) 
library(scater) 
library(scran) 
library(dplyr) 
library(scuttle) 
library(ggrepel) 
library(Seurat) 
library(ggplot2) 
library(ggbeeswarm) 
library(ggpubr) 
library(RColorBrewer) 
library(knitr) 

#BiocManager::install(c("limma", "edgeR"))
library(limma)
library(edgeR)
#BiocManager::install("miloR", lib = libpath)
library(miloR)
library(patchwork)


#Annotated object 
#all_first <- readRDS("/scratch/project_2012655/scRNAseq course/Downstream/all.rds")
all <- readRDS("/scratch/project_2012655/Aurora/all_donor_annotated_05.RDS")

#all$condition <- sapply(strsplit(all$orig.ident, "_"), function(x) paste(x[-1], collapse="_"))

#all.sub <- all[,!(all$condition=="vasc" | all$donor=="Alz37cl2")]
#all.sub <- all[,!(all$condition=="vasc")]


# to sce object
sce <- as.SingleCellExperiment(all)
colData(sce)$Sample <- as.character(all$orig.ident)        # näyte / replikaatti
#colData(sce)$condition <- as.character(all.sub$sampletype)     # condition (ctr/vasc/vasc_mg)


# jos sinulla on cell type -annotaatio:
#stopifnot(all(!is.null(all$sctype_pred)))

all_milo <- Milo(sce)


## stacked barplot (proportions)

# propTab <- table(colData(all_milo)$condition,  colData(all_milo)$sctype_pred)/rowSums(table(colData(all_milo)$condition,  colData(all_milo)$sctype_pred))
# propTab <- as.data.frame(propTab)
# colnames(propTab) <- c("Model","CellType","Proportion")
# 
# getPalette = colorRampPalette(brewer.pal(9, "Set1"))
# colourCount = length(unique(propTab$CellType))
# colVec_models <- setNames(getPalette(colourCount),
#                           sort(unique(propTab$CellType)))

# barplotStckd <- ggplot(propTab, aes(fill=CellType, x=Model, y=Proportion))+
#   geom_bar(position="fill", stat="identity")+
#   #facet_grid(. ~ unifiedSampleID, scales = "free_x", space = "free_x")+
#   theme_bw()+
#   theme(plot.title=element_text(size=18, face="bold", hjust=0.5),
#         legend.title = element_text(size = 11, face="bold", hjust=0.5),
#         legend.text  = element_text(size = 9),
#         legend.key.size = unit(0.5, "lines"),
#         axis.text.x = element_text(size=9, angle=90, vjust=0.5, hjust=1),
#         axis.text.y = element_text(size=9),
#         axis.title=element_text(size=11),
#         strip.text.x = element_text(size = 10)) +
#   guides(shape = guide_legend(override.aes = list(size = 5)),
#          fill = guide_legend(override.aes = list(size = 5), ncol=1))+
#   scale_fill_manual(name="Cell types",
#                     values = colVec_models)+
#   scale_y_continuous(labels = scales::percent, breaks=seq(0,1,0.2))+
#   xlab("Model")+
#   ylab("Cell type abundance [%]")
# 
dirQCplots <- "/scratch/project_2012655/Aurora/milor/"
# 
# pdf(paste0(dirQCplots,"barplotStckd_ctr_vasc_batch.pdf"))
# plot(barplotStckd)
# dev.off()
# 

## Barplot (numCells)

sampleCellsNum <- as.data.frame(table(colData(all_milo)$condition))
colnames(sampleCellsNum) <- c("Model", "Ncells")

barplotLines <- ggplot(sampleCellsNum, aes(x=Model, y=Ncells))+
  geom_bar(stat="identity", fill="black")+
  geom_text(aes(label=Ncells), vjust=-0.3, size=5)+
  theme_bw()+
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle=90, vjust=0.5, hjust=1),
        axis.title=element_text(size=9))+
  ylab("Number of cells")+xlab("Model")+
  scale_y_continuous(labels = scales::comma, breaks=seq(0,15000,2500), limits=c(0,15000))

barplotLines

# #SIMPLE CHECK
# tmp <- as.data.frame(table(sce$sctype_pred))
# colnames(tmp) <- c("Cluster", "Number of cells")
# kable(tmp)
# # UMAP-plotit



# plotUMAP <- plotReducedDim(all_milo, colour_by="sctype_pred", dimred = "UMAP", point_size=0.1)

plotUMAP2 <- plotReducedDim(all_milo, colour_by="condition",
                            dimred = "UMAP", point_size=0.3)+
  guides(colour = guide_legend(override.aes = list(size=10),
                               title="Clusters"))
plotUMAP2




plotUMAP3 <- DimPlot(
  all,
  group.by = "condition",
  reduction = "umap",
)+facet_wrap(~condition, nrow=2)

pdf(paste0(dirQCplots,"umapTestCondition_donor.pdf"))
plot(plotUMAP3)
dev.off()


plotUMAP4 <- DimPlot(
  all,
  group.by = "orig.ident",
  reduction = "umap",
)+facet_wrap(~orig.ident, nrow=3)

pdf(paste0(dirQCplots,"umapTestCondition_donor.pdf"))
plot(plotUMAP4)
dev.off()



#KNN graph and neighbourhoods

all_milo <- buildGraph(all_milo, k = 50, d = 30, reduced.dim = "HARMONY")
all_milo <- makeNhoods(all_milo, prop = 0.3, k = 50, d = 30, refined = TRUE, reduced_dims = "HARMONY")

#Neighbourhood plot
nhoodPlot <- plotNhoodSizeHist(all_milo)
nhoodPlot

all_milo <- countCells(all_milo, meta.data  = data.frame(colData(all_milo)), samples="Sample")

head(nhoodCounts(all_milo))

#Design 

colData(all_milo)$donor <- gsub("_.+","",colData(all_milo)$Sample)
design_df <- as.data.frame(colData(all_milo)[, c("Sample" ,"condition","donor")])

design_df <- distinct(design_df)

# Tee 1 rivi per sample
rownames(design_df) <- design_df$Sample



design_df <- design_df[colnames(nhoodCounts(all_milo)), , drop = FALSE]
design_df$condition <- factor(design_df$condition, levels=c("ctr","vasc","vasc_mg"))


message("Condition levels after ordering:")
print(table(design_df$condition))
if(length(levels(design_df$condition)) < 2) stop("Design sisältää vähemmän kuin 2 condition-tasoa; differential abundance ei voi toimia.")

## Option 1: Model all effects in contrast

all_milo <- calcNhoodDistance(all_milo, d = 30, reduced.dim = "HARMONY")

saveRDS(all_milo, "/scratch/project_2012655/Aurora/milor_Organoids_Annotated_donor_v7.RDS")



colData(all_milo)$donor <- factor(colData(all_milo)$donor, levels=c("Alz37cl2","Alz43cl2","PLCG1"))

contrast.3 <- c("conditionvasc - conditionctr",
                "conditionvasc_mg - conditionctr",
                "conditionvasc_mg - conditionvasc",
                "donorAlz43cl2",
                "donorPLCG1")

res_contrasts <- sapply(1:length(contrast.3), function(y){
  
  da_res <- testNhoods(
    all_milo,
    design = ~ 0 + condition + donor,
    design.df = design_df,
    model.contrasts = contrast.3[y],
    fdr.weighting="graph-overlap", norm.method="TMM"
  )
  
  da_res$contrast <- contrast.3[y]
  
  return(da_res)
  
}, simplify=F)

names(res_contrasts) <- contrast.3


# conditionvasc effects versus conditionvasc_mg effects
plot(res_contrasts[[1]]$logFC, res_contrasts[[2]]$logFC)
cor(res_contrasts[[1]]$logFC, res_contrasts[[2]]$logFC, method="pearson")
cor.test(res_contrasts[[1]]$logFC, res_contrasts[[2]]$logFC, method="pearson")$p.value

# DonorAlz43cl2 effects versus DonorPLCG1 effects
plot(res_contrasts[[4]]$logFC, res_contrasts[[5]]$logFC)
cor(res_contrasts[[4]]$logFC, res_contrasts[[5]]$logFC, method="pearson")
cor.test(res_contrasts[[4]]$logFC, res_contrasts[[5]]$logFC, method="pearson")$p.value

## comparison of conditionvasc effects vs DonorAlz43cl2 effects
plot(res_contrasts[[1]]$logFC, res_contrasts[[4]]$logFC)

## comparison of conditionvasc effects vs DonorAlz43cl2 effects
plot(res_contrasts[[1]]$logFC, (res_contrasts[[5]]$logFC))

da_res1 <- res_contrasts[[1]]
da_res2 <- res_contrasts[[2]]


ggplot(da_res1, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_res1, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) 

all_milo <- buildNhoodGraph(all_milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(all_milo, dimred = "UMAP", colour_by="condition", text_size = 3, point_size=0.5) + guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(all_milo, da_res1, layout="UMAP", alpha=0.1) 

umap_pl + nh_graph_pl +
  plot_layout(guides="collect")


all_milo <- buildNhoodGraph(all_milo)

## Find cell types

da_res1 <- annotateNhoods(all_milo, da_res1, coldata_col = "sctype_pred")
da_res2 <- annotateNhoods(all_milo, da_res2, coldata_col = "sctype_pred")

# Mixed-threshold
da_res1$sctype_pred <- ifelse(da_res1$sctype_pred_fraction < 0.7, "Mixed", da_res1$sctype_pred)
da_res2$sctype_pred <- ifelse(da_res2$sctype_pred_fraction < 0.7, "Mixed", da_res2$sctype_pred)

# (valinnainen) tarkista jakauma
ggplot(da_res1, aes(sctype_pred_fraction)) + geom_histogram(bins = 50)
ggplot(da_res2, aes(sctype_pred_fraction)) + geom_histogram(bins = 50)

# (valinnainen) beeswarm solutyypeittäin
plotDAbeeswarm(da_res1, group.by = "sctype_pred")
plotDAbeeswarm(da_res2, group.by = "sctype_pred")

library(ggplot2)
library(ggbeeswarm)

ggplot(da_res2, aes(x = sctype_pred, y = logFC)) +
  ggbeeswarm::geom_quasirandom(
    width = 0.2,
    size = 1.5,
    colour = "grey70"
  ) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  ) +
  labs(
    x = NULL,
    y = "logFC"
  )


# Lisää sama cell type -annotaatio kaikille kontrasteille
res_contrasts_annot <- lapply(res_contrasts, function(df) {
  df <- annotateNhoods(all_milo, df, coldata_col = "sctype_pred")
  df$sctype_pred <- ifelse(df$sctype_pred_fraction < 0.7, "Mixed", df$sctype_pred)
  df
})

allDA_annot <- do.call(rbind, res_contrasts_annot)
rownames(allDA_annot) <- NULL
allDA_annot$sctype_pred <- factor(allDA_annot$sctype_pred, levels=unique(allDA_annot$sctype_pred))

## Find groups
da_results <- groupNhoods(all_milo, da_res1, max.lfc.delta = 4)
head(da_results)


vecCor <- setNames(da_results$NhoodGroup, da_results$Nhood)


allDA_res <- sapply(res_contrasts, function(y){
  
  y$NhoodGroup <- vecCor[y$Nhood]
  return(y)
  
}, simplify=F)


allDA_res <- do.call("rbind", allDA_res); rownames(allDA_res) <- NULL

allDA_res$NhoodGroup <- factor(allDA_res$NhoodGroup, levels=unique(allDA_res$NhoodGroup))


##Violin plots

effectSize <- ggplot(allDA_res, aes(x = contrast, y = logFC, fill = contrast)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  facet_wrap(~ sctype_pred, scales = "free_y") +
  theme_bw()+
  theme(axis.text.x=element_blank())+
  ylab("logFC per contrast")+
  xlab("Contrast of interest per 'NhoodGroup'")+
  ggtitle("Analysing the effect size of the 3D model and the donor")


pdf(paste0(dirQCplots,"effectSize_contrasts.pdf"), width=12, height=6)
plot(effectSize)
dev.off()

df_plot <- allDA_annot %>%
  filter(!is.na(logFC), !is.na(contrast), !is.na(sctype_pred)) %>%
  group_by(contrast, sctype_pred) %>%
  mutate(n_grp = n()) %>%
  ungroup()

effectSize <- ggplot(df_plot, aes(x = contrast, y = logFC, fill = contrast)) +
  geom_violin(data = subset(df_plot, n_grp >= 2), trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  facet_wrap(~ sctype_pred, scales = "free_y", drop = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_blank())

pdf(paste0(dirQCplots,"effectSize_contrasts_annot.pdf"), width=12, height=6)
plot(effectSize)
dev.off()


sessionInfo()
library(uwot)

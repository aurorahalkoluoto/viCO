

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
library(limma)
library(edgeR)
library(miloR)
library(patchwork)







###############################
## 0. Pseudobulk per sample ## 
##############################

## Sample: Each donor in the three conditions (9 samples). Cell types are ignored
## The sample is treaten like a juice with all the fruits/cell types

## Check that at least all samples have 10 cells
stopifnot(!any(names(which(table(obj$orig.ident)<10))))

pseudo_obj <- AggregateExpression(
  obj,
  group.by = "orig.ident",  
  assays = "RNA",
  slot = "counts",  
  return.seurat=T
)


## Include information about batch, condition and donor

obj$orig.ident <- gsub("_","-",obj$orig.ident)

pseudo_obj$batch <- unname(obj$batch[match(pseudo_obj$orig.ident, obj$orig.ident)])
pseudo_obj$condition <- unname(obj$condition[match(pseudo_obj$orig.ident, obj$orig.ident)])
pseudo_obj$donor <- unname(obj$donor[match(pseudo_obj$orig.ident, obj$orig.ident)])



####################
## 1.PCA Fig S2C##
#################

## Find highly variable genes
pseudo_obj <- FindVariableFeatures(pseudo_obj)


matrixForPca <- prcomp(pseudo_obj@assays$RNA["data"][VariableFeatures(pseudo_obj),], scale = FALSE)
dfPca <- as.data.frame(matrixForPca$rotation)
stopifnot(all(stopifnot(rownames(pseudo_obj@meta.data)==rownames(dfPca))))

dfPca <- cbind(pseudo_obj@meta.data, dfPca)
dfPca$sampleId <- factor(dfPca$orig.ident, levels=unique(sort(dfPca$orig.ident)))

## Just to know how much the principle components of each PC explain (since we will plot just the first two)
explained_variance <- matrixForPca$sdev^2 / sum(matrixForPca$sdev^2)
explained_variance_non_scientific <- format(round(explained_variance*100,2), scientific = FALSE)
names(explained_variance_non_scientific) <- paste0("PC", 1:length(explained_variance_non_scientific))
print(explained_variance_non_scientific)


## Get colors for the samples (I think in your case they are 12)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(pseudo_obj$condition))
colVec <- setNames(getPalette(colourCount),
                   unique(pseudo_obj$condition))





##Plot just with two principal components: PC1-PC2
pseudoPCA <- ggplot(dfPca, aes(x=PC1, y=PC2, fill=condition, shape=donor, color=batch))+
  geom_point(size=5,stroke=1)+
  theme_bw()+
  scale_fill_manual(name="Condition", values=colVec)+
  scale_color_manual(name="Batch", values=c("black","blue"))+
  scale_shape_manual(name="Donor", values=c(21,22,23))+
  theme(plot.title=element_text(hjust=0.5, face="bold", size=14),
        legend.position="right")+
  ggtitle(paste0("Pseudobulked & HVG scaled"))+
  xlab(paste0("PC1, var.explained = ", explained_variance_non_scientific["PC1"], " %"))+
  ylab(paste0("PC2, var.explained = ", explained_variance_non_scientific["PC2"], " %"))+
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        color = "black"
      )
    ),
    color = guide_legend(
      override.aes = list(
        shape = 21,
        fill = "white", 
        stroke = 1
      )
    )
  )



dfPca <- dfPca %>%
  mutate(
    donor = recode(donor,
                   "Alz43cl2"  = "CL2",
                   "Alz37cl2" = "CL1",
                   "PLCG1"    = "CL3"),
    condition = recode(condition,
                       "ctr"     = "cCO",
                       "vasc"    = "vCO",
                       "vasc_mg" = "viCO")
  )

## --- Factor order (optional but recommended) ---
dfPca$donor <- factor(dfPca$donor, levels = c("CL1", "CL2", "CL3"))
dfPca$condition <- factor(dfPca$condition, levels = c("cCO", "vCO", "viCO"))

## --- Colors for condition (fill) with alpha ---
colVec <- c(
  "cCO"  = scales::alpha("#0072B2", 0.5),
  "vCO"  = scales::alpha("#D55E00", 0.5),
  "viCO" = scales::alpha("#7B3294", 0.5)
)

## --- Plot: PC1 vs PC2 ---
pseudoPCA <- ggplot(dfPca, aes(x = PC1, y = PC2, fill = condition, shape = donor)) +
  geom_point(size = 5, stroke = 1) +
  theme_bw() +
  scale_fill_manual(name = "Condition", values = colVec) +
  scale_shape_manual(name = "Donor", values = c(21, 22, 23)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  ) +
  ggtitle("Pseudobulked & HVG scaled") +
  xlab(paste0("PC1, var.explained = ", explained_variance_non_scientific["PC1"], " %")) +
  ylab(paste0("PC2, var.explained = ", explained_variance_non_scientific["PC2"], " %")) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, color = "black")
    ),
    color = guide_legend(
      override.aes = list(shape = 21, fill = "white", stroke = 1)
    )
  )





## grid with PCA correlation

pc_pairs <- combn(paste0("PC", 1:6), 2, simplify = FALSE)

library(dplyr)
library(purrr)

plot_df <- map_dfr(pc_pairs, function(pcs) {
  dfPca %>%
    transmute(
      PCx = .data[[pcs[1]]],
      PCy = .data[[pcs[2]]],
      PCx_name = pcs[1],
      PCy_name = pcs[2]
    )
})

plot_df$sampleId <- gsub("\\..+","",rownames(plot_df)) 
plot_df$condition <- dfPca[match(plot_df$sampleId, dfPca$sampleId),]$condition
plot_df$donor <- dfPca[match(plot_df$sampleId, dfPca$sampleId),]$donor


gridpca_first6 <- ggplot(plot_df, aes(PCx, PCy, col=condition, shape=donor)) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_grid(
    PCy_name ~ PCx_name,
    scales = "free",
    drop = TRUE
  ) +
  theme_bw()+
  theme(axis.text.x=element_text(hjust=1, vjust=0.5, size=8, angle=90),
        axis.text.y=element_text(size=8))



pc_pairs <- combn(paste0("PC", 1:6), 2, simplify = FALSE)



library(dplyr)
library(purrr)
library(ggplot2)
library(scales)

plot_df <- map_dfr(pc_pairs, function(pcs) {
  dfPca %>%
    transmute(
      PCx = .data[[pcs[1]]],
      PCy = .data[[pcs[2]]],
      PCx_name = pcs[1],
      PCy_name = pcs[2]
    )
})

plot_df$sampleId <- gsub("\\..+","", rownames(plot_df))
plot_df$condition <- dfPca[match(plot_df$sampleId, dfPca$sampleId), ]$condition
plot_df$donor     <- dfPca[match(plot_df$sampleId, dfPca$sampleId), ]$donor

# Factor order (suositeltava)
plot_df$condition <- factor(plot_df$condition, levels = c("cCO","vCO","viCO"))
plot_df$donor <- factor(plot_df$donor, levels = c("CL1","CL2","CL3"))

# Samat värit kuin aiemmin (mutta color-estetiikkaan ilman alpha-haalistusta)
colVec_line <- c(
  "cCO"  = "#0072B2",
  "vCO"  = "#D55E00",
  "viCO" = "#7B3294"
)

gridpca_first6 <- ggplot(plot_df, aes(PCx, PCy, col = condition, shape = donor)) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_grid(
    PCy_name ~ PCx_name,
    scales = "free",
    drop = TRUE
  ) +
  theme_bw() +
  scale_color_manual(name = "Condition", values = colVec_line) +
  theme(
    axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 8, angle = 90),
    axis.text.y = element_text(size = 8)
  )



##################################
### 2. GSVA for the pathways   ##
### Fig 3D,S3A (all) & S4B(RG) ##
#################################

# Continue with pseudobulked sample with all cells or RG (From stem 0)


library(GSVA)
library(ggbeeswarm)
#library(gghighlight)
library(ggpubr)
library(tidyverse)
library(GSEABase)


# Load the GMT file
gmt_file_c2 <- "/scratch/project_2012655/OT/msigdb/c2.all.v2026.1.Hs.symbols.gmt" 
gmt_file_h <- "/scratch/project_2012655/OT/msigdb/h.all.v2026.1.Hs.symbols.gmt" 
gmt_file_c5 <- "/scratch/project_2012655/OT/msigdb/c5.go.bp.v2026.1.Hs.symbols.gmt" 
gene_sets_c2 <- getGmt(gmt_file_c2) 
gene_sets_c5 <- getGmt(gmt_file_c5)
gene_sets_h <- getGmt(gmt_file_h)

relevant_pathways <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1",
  "KEGG_CELL_CYCLE",
  "KEGG_DNA_REPLICATION",
  "HALLMARK_DNA_REPAIR"
)


c2_list <- sapply(relevant_pathways[relevant_pathways %in% names(gene_sets_c2)],
                  function(y){
                    
                    geneIds(gene_sets_c2[[y]])
                    
                  }, simplify=F)

c5_list <- sapply(relevant_pathways[relevant_pathways %in% names(gene_sets_c5)],
                  function(y){
                    
                    geneIds(gene_sets_c5[[y]])
                    
                  }, simplify=F)

h_list <- sapply(relevant_pathways[relevant_pathways %in% names(gene_sets_h)],
                 function(y){
                   
                   geneIds(gene_sets_h[[y]])
                   
                 }, simplify=F)


pathways <- c(c2_list, h_list, c5_list)
pathways_filt <- sapply(pathways, function(x) x[x %in% rownames(pseudo_obj@assays$RNA["data"])])


gsvaPar <- gsvaParam(pseudo_obj@assays$RNA["data"], pathways_filt)
gsva.es <- gsva(gsvaPar, verbose=FALSE)
dim(gsva.es)
fullMatrix_results <- t(gsva.es[1:dim(gsva.es)[1],1:dim(gsva.es)[2]])
fullMatrix_results <- as.data.frame(fullMatrix_results)
fullMatrix_results$sampleId <- rownames(fullMatrix_results)
rownames(fullMatrix_results) <- NULL



gsva_long <- as.data.frame(fullMatrix_results %>% pivot_longer(-c("sampleId"), names_to="Pathways", values_to="GSVA_scores"))

pathway_groups <- list(
  "stem" = c(
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1",
    "KEGG_CELL_CYCLE",
    "KEGG_DNA_REPLICATION",
    "HALLMARK_DNA_REPAIR"))

pathway_to_group <- stack(pathway_groups)
colnames(pathway_to_group) <- c("Pathways", "PathwayGroup")

gsva_long$PathwayGroup <- pathway_to_group$PathwayGroup[
  match(gsva_long$Pathways, pathway_to_group$Pathways)
]

pathway_order <- unlist(pathway_groups, use.names = FALSE)
gsva_long$Pathways <- factor(gsva_long$Pathways, levels = pathway_order)
gsva_long$PathwayGroup <- factor(gsva_long$PathwayGroup, levels = names(pathway_groups))


obj$orig.ident <- gsub("_","-",obj$orig.ident)

gsva_long$batch <- unname(obj$batch[match(gsva_long$sampleId, obj$orig.ident)])
gsva_long$condition <- unname(obj$condition[match(gsva_long$sampleId, obj$orig.ident)])
gsva_long$donor <- unname(obj$donor[match(gsva_long$sampleId, obj$orig.ident)])

gsva_long$condition <- factor(gsva_long$condition, levels=c("ctr","vasc","vasc_mg"))

gsva_long$condition <- dplyr::recode(
  gsva_long$condition,
  "ctr"     = "cCO",
  "vasc"    = "vCO",
  "vasc_mg" = "viCO"
)

gsva_long$donor <- dplyr::recode(
  gsva_long$donor,
  "Alz37cl2" = "CL1",
  "Alz43cl2" = "CL2",
  "PLCG1"    = "CL3"
)

## Set plotting order (optional but recommended)
gsva_long$condition <- factor(gsva_long$condition, levels = c("cCO","vCO","viCO"))
gsva_long$donor     <- factor(gsva_long$donor, levels = c("CL1","CL2","CL3"))
gsva_long$sampleId_pretty <- paste(gsva_long$condition)

allrelevPaths_gsva <- ggplot(gsva_long, aes(y=Pathways, x=sampleId_pretty, fill=GSVA_scores))+
  geom_tile(col="black")+
  facet_grid(PathwayGroup ~ donor, scales="free", space="free")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8),
        axis.text.y=element_text(size=4),
        strip.text=element_text(size=7),
        strip.text.y.right = element_blank(),
        legend.position="top")+
  xlab("Activation Scores")+
  scale_fill_gradient2(midpoint=0, low="blue", mid="white",
                       high="red",
                       guide = guide_colorbar(
                         title.position = "top",
                         title.hjust = 0.5,
                         barwidth = 10,
                         barheight = 0.8
                       ))+
  ylab("")

pdf(file=paste0(dirQCplots,"gsvaScores_allpathways_RGa.pdf"), width=12, height=4)
plot(allrelevPaths_gsva)
dev.off()



allrelevPaths_gsva2 <- ggplot(gsva_long, aes(x = condition, y = GSVA_scores, colour = donor)) +
  facet_wrap(~ Pathways) +
  ggtitle("Activation scores per pathway")+
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.1,   # jitter within condition
      dodge.width  = 0.5     # separation between conditions
    ),
    alpha = 0.8
  ) +
  theme_bw()+
  theme(legend.position="top",
        strip.text=element_text(size=4.5))
  # scale_color_manual(values=c("red","blue"))


pdf(file=paste0(dirQCplots,"gsvaScores_point_allpathways_sabita_RGb.pdf"), width=11, height=8)
plot(allrelevPaths_gsva2)
dev.off()


#################
# Fig 3C #######
################


#VENN + variance

library(Seurat)
library(dplyr)
library(stringr)
library(VennDiagram)
library(ggplot2)


stopifnot(all(c("condition", "orig.ident") %in% colnames(all@meta.data)))

RG <- subset(all, subset = sctype_pred == "Radial glia")
RG <- subset(all, subset = sctype_pred == "Neuroblast")



### iterate over "organoid model, aka condition"
organoid_models <- unique(RG$condition)

res_model <- sapply(organoid_models, function(z){
  
  tmp_subset <- subset(RG, condition==z)
  
  vecNames <- setNames(c("CL1","CL2","CL3"),
                       c("Alz37cl2", "Alz43cl2", "PLCG1"))
  
  tmp_subset$correctName <- unname(vecNames[tmp_subset$donor])
  lines_to_compare <- unique(tmp_subset$correctName)

  tmp_subset$cl1_specific <- ifelse(tmp_subset$correctName=="CL1","CL1","Others")
  tmp_subset$cl2_specific <- ifelse(tmp_subset$correctName=="CL2","CL2","Others")
  tmp_subset$cl3_specific <- ifelse(tmp_subset$correctName=="CL3","CL3","Others")

  
  run_deg <- function(obj, g1, g2, fdr = 0.05, lfc = 0.25,
                      test.use = "wilcox", min.pct = 0) {
    
    res <- FindMarkers(
      object = obj,
      ident.1 = g1,
      ident.2 = g2,
      test.use = test.use,
      logfc.threshold = 0,   
      min.pct = min.pct
    )
    
    res$gene <- rownames(res)
    res$test <- g1
    res$organoidModel <- z
      
    return(res)

  }
  
  Idents(tmp_subset) <- "cl1_specific"
  deg_cl1 <- run_deg(tmp_subset, "CL1", "Others")
  
  Idents(tmp_subset) <- "cl2_specific"
  deg_cl2 <- run_deg(tmp_subset, "CL2", "Others")
  
  Idents(tmp_subset) <- "cl3_specific"
  deg_cl3 <- run_deg(tmp_subset, "CL3", "Others")
  

  deg_model <- do.call("rbind", list(deg_cl1, deg_cl2, deg_cl3))
  
  return(deg_model)
  
}, simplify=F)


## get dysregulated all together

get_genes_dysregulated <- function(line="CL1", model="ctr"){
  
  lfc_threshold <- 0.25
  fdr = 0.05
  tmp <- res_model[[model]]
  tmp <- subset(tmp, test == line)
  signifGenes <- subset(tmp, p_val_adj<fdr & abs(avg_log2FC)>lfc_threshold)$gene

  return(signifGenes)
  
}


venn_plot_function <- function(venn_sets_named, cl1="#08306B",cl2="#2171B5",cl3="#6BAED6", nameFile="test.pdf"){
  
  venn_plot <- venn.diagram(
    x = venn_sets_named,
    filename = NULL,
    fill = c(cl1, cl2, cl3),
    alpha = 0.6,
    col = "black",
    lwd = 1.2,
    
    fontfamily = "sans",
    cat.fontfamily = "sans",
    
    cex = 1.6,
    cat.cex = 1.5,
    
    
    cat.dist = c(0.07, 0.07, 0.07),
    
    margin = 0.1
  )
  
  grid.newpage()
  grid.draw(venn_plot)
  
  # --- 7) (valinnainen) tallenna tulokset
  dir_out <- "/scratch/project_2012655/Aurora/plots/RG"
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file.path(dir_out, nameFile), venn_plot, width = 7, height = 6)
}
  
  

## cCO (RG) -dysregulated

venn_sets_named <- list(
  CL1=get_genes_dysregulated(line="CL1", model="ctr"),
  CL2=get_genes_dysregulated(line="CL2", model="ctr"),
  CL3=get_genes_dysregulated(line="CL3", model="ctr")
)

venn_sets_named <- lapply(venn_sets_named, function(x) {
  unique(na.omit(as.character(x)))
})

venn_plot_function(venn_sets_named, cl1="#08306B",cl2="#2171B5",cl3="#6BAED6", nameFile="Venn_CTR_celllineRG_DEGs.pdf")
  

## vCO (RG) -dysregulated
venn_sets_named <- list(
  CL1=get_genes_dysregulated(line="CL1", model="vasc"),
  CL2=get_genes_dysregulated(line="CL2", model="vasc"),
  CL3=get_genes_dysregulated(line="CL3", model="vasc")
)

venn_sets_named <- lapply(venn_sets_named, function(x) {
  unique(na.omit(as.character(x)))
})

venn_plot_function(venn_sets_named, cl1="#fdd7b5",cl2="#e89166",cl3="#b17d68", nameFile="Venn_vCO_celllineRG_DEGs.pdf")


## viCO (RG) -dysregulated

venn_sets_named <- list(
  CL1=get_genes_dysregulated(line="CL1", model="vasc_mg"),
  CL2=get_genes_dysregulated(line="CL2", model="vasc_mg"),
  CL3=get_genes_dysregulated(line="CL3", model="vasc_mg")
)

venn_sets_named <- lapply(venn_sets_named, function(x) {
  unique(na.omit(as.character(x)))
})

venn_plot_function(venn_sets_named, cl1="#8b65b0",cl2="#a696c7",cl3="#d7d8e9", nameFile="Venn_viCO_celllineRG_DEGs.pdf")


venn_plot_function_direction <- function(venn_sets_named, cl1="#08306B",cl2="#2171B5",cl3="#6BAED6", nameFile="test.pdf"){
  
  venn_plot <- venn.diagram(
    x = venn_sets_named,
    filename = NULL,
    fill = c(cl1, cl2, cl3),
    alpha = 0.6,
    col = "black",
    lwd = 1.2,
    
    fontfamily = "sans",
    cat.fontfamily = "sans",
    
    cex = 1.6,
    cat.cex = 0
    
  )
  
  grid.newpage()
  grid.draw(venn_plot)
  grid.text("CL1", x = 0.2, y = 0.75,
            gp = gpar(fontsize = 16, fontfamily = "sans", fontface = "bold"))
  
  grid.text("CL2", x = 0.5, y = 0.75,
            gp = gpar(fontsize = 16, fontfamily = "sans", fontface = "bold"))
  
  grid.text("CL3", x = 0.8, y = 0.75,
            gp = gpar(fontsize = 16, fontfamily = "sans", fontface = "bold"))
  dir_out <- "/scratch/project_2012655/Aurora/plots/RG"
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file.path(dir_out, nameFile), venn_plot, width = 7, height = 6)
}



#################################################################
## make the plots with upregulated and downregulated separated ##
#################################################################

get_genes_upregulated <- function(line="CL1", model="ctr"){
  lfc_threshold <- 0.25
  fdr = 0.05
  tmp <- res_model[[model]]
  tmp <- subset(tmp, test==line)
  signifGenes <- subset(tmp, p_val_adj<fdr & avg_log2FC>lfc_threshold)$gene
  return(signifGenes)
}

get_genes_downregulated <- function(line="CL1", model="ctr"){
  lfc_threshold <- 0.25
  fdr = 0.05
  tmp <- res_model[[model]]
  tmp <- subset(tmp, test==line)
  signifGenes <- subset(tmp, p_val_adj<fdr & avg_log2FC < -lfc_threshold)$gene
  return(signifGenes)
}


## cCO (RG) -upregulated
venn_sets_named_up <- list(
  CL1=get_genes_upregulated(line="CL1", model="ctr"),
  CL2=get_genes_upregulated(line="CL2", model="ctr"),
  CL3=get_genes_upregulated(line="CL3", model="ctr")
)

## cCO (RG) -downregulated
venn_sets_named_down <- list(
  CL1=get_genes_downregulated(line="CL1", model="ctr"),
  CL2=get_genes_downregulated(line="CL2", model="ctr"),
  CL3=get_genes_downregulated(line="CL3", model="ctr")
)


venn_sets_named_up <- lapply(venn_sets_named_up, function(x) {
  unique(na.omit(as.character(x)))
})

venn_sets_named_down <- lapply(venn_sets_named_down, function(x) {
  unique(na.omit(as.character(x)))
})


venn_plot_function_direction(venn_sets_named_up, cl1="#08306B",cl2="#2171B5",cl3="#6BAED6", nameFile="Venn_CTR_celllineRG_DEGs_upregulated.pdf")
venn_plot_function(venn_sets_named_down, cl1="#08306B",cl2="#2171B5",cl3="#6BAED6", nameFile="Venn_CTR_celllineRG_DEGs_downregulated.pdf")



## vCO (RG) -upregulated
venn_sets_named_up <- list(
  CL1=get_genes_upregulated(line="CL1", model="vasc"),
  CL2=get_genes_upregulated(line="CL2", model="vasc"),
  CL3=get_genes_upregulated(line="CL3", model="vasc")
)

## vCO (RG) -downregulated
venn_sets_named_down <- list(
  CL1=get_genes_downregulated(line="CL1", model="vasc"),
  CL2=get_genes_downregulated(line="CL2", model="vasc"),
  CL3=get_genes_downregulated(line="CL3", model="vasc")
)


venn_sets_named_up <- lapply(venn_sets_named_up, function(x) {
  unique(na.omit(as.character(x)))
})

venn_sets_named_down <- lapply(venn_sets_named_down, function(x) {
  unique(na.omit(as.character(x)))
})


venn_plot_function(venn_sets_named_up, cl1="#fdd7b5",cl2="#e89166",cl3="#b17d68", nameFile="Venn_vCO_celllineRG_DEGs_upregulated.pdf")
venn_plot_function(venn_sets_named_down, cl1="#fdd7b5",cl2="#e89166",cl3="#b17d68", nameFile="Venn_vCO_celllineRG_DEGs_downregulated.pdf")



# viCO (RG) -upregulated
venn_sets_named_up <- list(
  CL1=get_genes_upregulated(line="CL1", model="vasc_mg"),
  CL2=get_genes_upregulated(line="CL2", model="vasc_mg"),
  CL3=get_genes_upregulated(line="CL3", model="vasc_mg")
)

## viCO (RG) -downregulated
venn_sets_named_down <- list(
  CL1=get_genes_downregulated(line="CL1", model="vasc_mg"),
  CL2=get_genes_downregulated(line="CL2", model="vasc_mg"),
  CL3=get_genes_downregulated(line="CL3", model="vasc_mg")
)


venn_sets_named_up <- lapply(venn_sets_named_up, function(x) {
  unique(na.omit(as.character(x)))
})

venn_sets_named_down <- lapply(venn_sets_named_down, function(x) {
  unique(na.omit(as.character(x)))
})


venn_plot_function_direction(venn_sets_named_up, cl1="#8b65b0",cl2="#a696c7",cl3="#d7d8e9", nameFile="Venn_viCO_celllineRG_DEGs_upregulated.pdf")
venn_plot_function(venn_sets_named_down, cl1="#8b65b0",cl2="#a696c7",cl3="#d7d8e9",nameFile="Venn_viCO_celllineRG_DEGs_downregulated.pdf")



# GO 

library(clusterProfiler)
library(org.Hs.eg.db)


run_GO <- function(genes){
  
  genes <- unique(na.omit(as.character(genes)))
  
  if(length(genes) < 10){
    message("Too few genes, skipping GO")
    return(NULL)
  }
  
  gene_df <- tryCatch({
    bitr(genes,
         fromType = "SYMBOL",
         toType = "ENTREZID",
         OrgDb = org.Hs.eg.db)
  }, error = function(e) return(NULL))
  
  if(is.null(gene_df) || nrow(gene_df) < 10){
    message("Gene conversion failed or too few mapped genes")
    return(NULL)
  }
  
  ego <- tryCatch({
    enrichGO(
      gene          = gene_df$ENTREZID,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      readable      = TRUE
    )
  }, error = function(e) return(NULL))
  
  return(ego)
}




conditions <- c("ctr", "vasc", "vasc_mg")
lines <- c("CL1", "CL2", "CL3")

GO_results <- list()

for(cond in conditions){
  
  GO_results[[cond]] <- list()
  
  for(cl in lines){
    
    cat("Running GO:", cond, cl, "\n")
    
    genes_up <- get_genes_upregulated(cl, cond)
    genes_down <- get_genes_downregulated(cl, cond)
    
    GO_results[[cond]][[cl]] <- list(
      up = run_GO(genes_up),
      down = run_GO(genes_down)
    )
  }
}


dir_go <- "/scratch/project_2012655/Aurora/plots/RG/GO"
dir.create(dir_go, recursive = TRUE, showWarnings = FALSE)

for(cond in names(GO_results)){
  for(cl in names(GO_results[[cond]])){
    for(direction in c("up","down")){
      
      ego <- GO_results[[cond]][[cl]][[direction]]
      
      if(!is.null(ego)){
        fname <- paste0("GO_", cond, "_", cl, "_", direction, ".csv")
        write.csv(as.data.frame(ego),
                  file = file.path(dir_go, fname),
                  row.names = FALSE)
      }
    }
  }
}


condition_map <- c(
  ctr = "cCO",
  vasc = "vCO",
  vasc_mg = "viCO"
)

dir_plot <- "/scratch/project_2012655/Aurora/plots/RG/GO"
dir.create(dir_plot, recursive = TRUE, showWarnings = FALSE)

for(cond in names(GO_results)){
  for(cl in names(GO_results[[cond]])){
    for(direction in c("up","down")){
    
      ego <- GO_results[[cond]][[cl]][[direction]]
      
      if(!is.null(ego) && nrow(as.data.frame(ego)) > 0){
        
        pdf(file.path(dir_plot,
                      paste0("GO_", condition_map[cond], "_", cl, "_", direction, ".pdf")),
            width = 7, height = 6)
        
        tryCatch({
          print(dotplot(ego, showCategory = 15) +
                  ggtitle(paste(condition_map[cond], cl, direction)))
        }, error = function(e){
          message("Plot failed: ", cond, cl, direction)
        })
        
        dev.off()
      }
    }
  }
}




















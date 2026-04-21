

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
library(limma)
library(edgeR)
library(miloR)
library(patchwork)

obj <- readRDS("/scratch/project_2012655/Aurora/all_donor_annotated_05.RDS")


obj2 <- readRDS("/scratch/project_2012655/Aurora/all_donor_annotated_05_RGsub.RDS")
stopifnot(all(colnames(obj)==colnames(obj2)))
obj$RG_commitment <- obj2$RG_commitment
rm(obj2)
gc()


ctypes_to_discard <- table(obj$orig.ident, obj$sctype_pred)>=10
ctypes_to_discard <- names(which(colSums(ctypes_to_discard)<dim(ctypes_to_discard)[1]))

obj_sub <- obj[,!obj$sctype_pred %in% ctypes_to_discard]

stopifnot(all(table(obj_sub$orig.ident, obj_sub$sctype_pred)>=10))

## Pseudobulk by sample

pseudo_obj_sub <- AggregateExpression(
  obj_sub,
  group.by = c("orig.ident","sctype_pred"),  # Replace with your grouping variable (e.g., sample, condition)
  assays = "RNA",
  slot = "counts",  # Use log-normalized expression instead of raw counts
  return.seurat=T
)

## Find highly variable genes
pseudo_obj_sub <- FindVariableFeatures(pseudo_obj_sub)


## Include information about batch, condition and donor

pseudo_obj_sub$old.ident <- sapply(strsplit(pseudo_obj_sub$orig.ident,"_"), function(x) x[1])
obj$orig.ident <- gsub("_","-", obj$orig.ident)

pseudo_obj_sub$batch <- unname(obj$batch[match(pseudo_obj_sub$old.ident, obj$orig.ident)])
pseudo_obj_sub$condition <- unname(obj$condition[match(pseudo_obj_sub$old.ident, obj$orig.ident)])
pseudo_obj_sub$donor <- unname(obj$donor[match(pseudo_obj_sub$old.ident, obj$orig.ident)])


#######_________________________________________________________#######

#####     1.     example test for Radial glia (DEGs+GO)
######________________________________________________________#######


pseudo_obj_sub_RG <- subset(pseudo_obj_sub, sctype_pred=="Radial glia")
dirQCplots <- "/scratch/project_2012655/Aurora/plots/RG/"


#vascular <- c("Endothelial", "Smooth muscle", "Fibroblast")
#pseudo_obj_sub_RG <- subset(
 # pseudo_obj_sub,
#  sctype_pred %in% vascular)

library(SummarizedExperiment)
library(edgeR)
library(limma)
library(sva)

#### differential expression

se <- SummarizedExperiment(assays=list(counts=as.matrix(pseudo_obj_sub_RG@assays$RNA["counts"])),
                           colData=pseudo_obj_sub_RG@meta.data)

dge <- DGEList(counts=assays(se)$counts, genes=as.data.frame(mcols(se)))
dim(dge)

assays(se)$logCPM <- edgeR::cpm(dge, log=TRUE, prior.count=0.5)

## Check the distribution

avgexp <- rowMeans(assays(se)$logCPM)

mask <- avgexp > -2
dim(se)
se.filt <- se[mask, ]
dim(se.filt)

dim(dge)
dge.filt <- dge[mask, ]
dim(dge.filt)

## Copy without normalisation of dge.filt counts to use later for voomQualityWeights
dge.filt.noTMM <- dge.filt

## Renormalise again
assays(se.filt)$logCPM <- edgeR::cpm(dge.filt, normalized.lib.sizes=TRUE,
                                     log=TRUE, prior.count=0.5)

colData(se.filt)$condition <- factor(colData(se.filt)$condition, levels=c("ctr","vasc","vasc_mg"))
colData(se.filt)$batch <- factor(colData(se.filt)$batch, levels=c("B1","B2"))
colData(se.filt)$donor <- factor(colData(se.filt)$donor, levels=c("Alz37cl2","Alz43cl2","PLCG1"))

se.filt$condition <- relevel(se.filt$condition, ref="ctr")

mod <- model.matrix(~ condition + donor, colData(se.filt))
mod0 <- model.matrix(~ 1, colData(se.filt))

IQRs <- apply(assays(se.filt)$logCPM, 1, IQR)
sv <- sva(assays(se.filt)$logCPM[IQRs > quantile(IQRs, prob=0.9), ],
          mod=mod, mod0=mod0)

cnames <- c(colnames(mod), paste0("SV", 1:sv$n))
mod <- cbind(mod, sv$sv)
colnames(mod) <- cnames
head(mod, n=3)

dge.filt <- estimateDisp(dge.filt, mod, robust=TRUE)
prdf <- cut(dge.filt$prior.df, breaks=c(0, 1, 2, 3, 4, 5))
table(prdf, useNA="always")

v <- voomWithQualityWeights(dge.filt.noTMM, mod, normalize="quantile")

fit <- lmFit(v, mod)
fit <- eBayes(fit, robust=TRUE)

res <- decideTests(fit, p.value=0.1)
summary(res)
#vascular <- c("Endothelial", "Smooth muscle", "Fibroblast")
#pseudo_obj_sub_RG <- subset(
 # pseudo_obj_sub,
#  sctype_pred %in% vascular)

library(SummarizedExperiment)
library(edgeR)
library(limma)
library(sva)

#### differential expression

se <- SummarizedExperiment(assays=list(counts=as.matrix(pseudo_obj_sub_RG@assays$RNA["counts"])),
                           colData=pseudo_obj_sub_RG@meta.data)

dge <- DGEList(counts=assays(se)$counts, genes=as.data.frame(mcols(se)))
dim(dge)

assays(se)$logCPM <- edgeR::cpm(dge, log=TRUE, prior.count=0.5)

## Check the distribution

avgexp <- rowMeans(assays(se)$logCPM)

mask <- avgexp > -2
dim(se)
se.filt <- se[mask, ]
dim(se.filt)

dim(dge)
dge.filt <- dge[mask, ]
dim(dge.filt)

## Copy without normalisation of dge.filt counts to use later for voomQualityWeights
dge.filt.noTMM <- dge.filt

## Renormalise again
assays(se.filt)$logCPM <- edgeR::cpm(dge.filt, normalized.lib.sizes=TRUE,
                                     log=TRUE, prior.count=0.5)

colData(se.filt)$condition <- factor(colData(se.filt)$condition, levels=c("ctr","vasc","vasc_mg"))
colData(se.filt)$batch <- factor(colData(se.filt)$batch, levels=c("B1","B2"))
colData(se.filt)$donor <- factor(colData(se.filt)$donor, levels=c("Alz37cl2","Alz43cl2","PLCG1"))

se.filt$condition <- relevel(se.filt$condition, ref="ctr")

mod <- model.matrix(~ condition + donor, colData(se.filt))
mod0 <- model.matrix(~ 1, colData(se.filt))

IQRs <- apply(assays(se.filt)$logCPM, 1, IQR)
sv <- sva(assays(se.filt)$logCPM[IQRs > quantile(IQRs, prob=0.9), ],
          mod=mod, mod0=mod0)

cnames <- c(colnames(mod), paste0("SV", 1:sv$n))
mod <- cbind(mod, sv$sv)
colnames(mod) <- cnames
head(mod, n=3)

dge.filt <- estimateDisp(dge.filt, mod, robust=TRUE)
prdf <- cut(dge.filt$prior.df, breaks=c(0, 1, 2, 3, 4, 5))
table(prdf, useNA="always")

v <- voomWithQualityWeights(dge.filt.noTMM, mod, normalize="quantile")

fit <- lmFit(v, mod)
fit <- eBayes(fit, robust=TRUE)

res <- decideTests(fit, p.value=0.1)
summary(res)

#        (Intercept) conditionvasc conditionvasc_mg donorAlz43cl2 donorPLCG1   SV1   SV2
# Down          1072             7               51          4109       3618  1294    14
# NotSig        3877         19277            19217         11254      11995 16731 19201
# Up           14343             8               24          3929       3679  1267    77

df_summary <- as.data.frame(summary(res))
colnames(df_summary) <- c("Direction","Covariate","numDE")

### Plot per covariate

resFit_conditionvasc <- topTable(fit, coef="conditionvasc", n=Inf)
resFit_conditionvasc_mg <- topTable(fit, coef="conditionvasc_mg", n=Inf)

# upReg_list <- rownames(resFit_disease_statusPatient[which(resFit_disease_statusPatient$adj.P.Val<0.1 & resFit_disease_statusPatient$logFC>0),])
# downReg_list <- rownames(resFit_disease_statusPatient[which(resFit_disease_statusPatient$adj.P.Val<0.1 & resFit_disease_statusPatient$logFC<0),])


### volcanoPlot

library(ggrepel)

resFit_conditionvasc$signif <- ifelse(resFit_conditionvasc$adj.P.Val<0.1, TRUE, FALSE)
resFit_conditionvasc$gene <- rownames(resFit_conditionvasc)

volcanoPlot <- ggplot(data=resFit_conditionvasc, aes(x=logFC, y=-log10(adj.P.Val), col=signif))+
  geom_point(size=1, pch=15)+
  theme_bw()+
  scale_color_manual(name="DE genes",values=c("TRUE"="red","FALSE"="grey80"))+
  geom_hline(yintercept=-log10(0.1), linetype="dashed", col="grey30")+
  geom_vline(xintercept=c(-log2(1.5),-log2(2),-log2(3),log2(1.5), log2(2), log2(3), log2(4)), linetype="dashed", col="grey30")+
  geom_text_repel(data=subset(resFit_conditionvasc, adj.P.Val<0.1), aes(label=gene), col="black")+
  xlab(expression(log[2]*"Fold Change")) +
  ylab(expression(-log[10]("adj. p-value"))) +
  ggtitle("vCO coeff., with SVA")+
  annotate("text", x = -0.13, y = -log10(0.12), label = "10% FDR", size = 4, hjust=0)+
  annotate("text", x = -log2(1.5)+0.1, y = 0, label = "<1/1.5 FC", size = 3, hjust=0, vjust=0.5, angle=90)+
  annotate("text", x = -log2(2)+0.1, y = 0, label = "<1/2 FC", size = 3, hjust=0, vjust=0.5,  angle=90)+
  annotate("text", x = -log2(3)+0.1, y = 0, label = "<1/3 FC", size = 3, hjust=0, vjust=0.5, angle=90)+
  annotate("text", x = log2(1.5)+0.1, y = 0, label = ">1.5 FC", size = 3, hjust=0, vjust=0.5, angle=90)+
  annotate("text", x = log2(2)+0.1, y = 0, label = ">2 FC", size = 3, hjust=0, vjust=0.5,  angle=90)+
  annotate("text", x = log2(3)+0.1, y = 0, label = ">3 FC", size = 3, hjust=0, vjust=0.5, angle=90)+
  theme(plot.title = element_text(hjust=0.5, size=16),
        legend.key.size = unit(1.5, "lines"),
        axis.text=element_text(size=14)
  )+
  guides(
    color = guide_legend(override.aes = list(size = 6)),
  )


pdf(file=paste0(dirQCplots,"volcanoPlot_vCO_RG.pdf"), width=8, height=7)
plot(volcanoPlot)
dev.off()


resFit_conditionvasc_mg$signif <- ifelse(resFit_conditionvasc_mg$adj.P.Val<0.1, TRUE, FALSE)
resFit_conditionvasc_mg$gene <- rownames(resFit_conditionvasc_mg)

volcanoPlot <- ggplot(data=resFit_conditionvasc_mg, aes(x=logFC, y=-log10(adj.P.Val), col=signif))+
  geom_point(size=1, pch=15)+
  theme_bw()+
  scale_color_manual(name="DE genes",values=c("TRUE"="red","FALSE"="grey80"))+
  geom_hline(yintercept=-log10(0.1), linetype="dashed", col="grey30")+
  geom_vline(xintercept=c(-log2(1.5),-log2(2),-log2(3),log2(1.5), log2(2), log2(3), log2(4)), linetype="dashed", col="grey30")+
  geom_text_repel(data=subset(resFit_conditionvasc_mg, adj.P.Val<0.1), aes(label=gene), col="black")+
  xlab(expression(log[2]*"Fold Change")) +
  ylab(expression(-log[10]("adj. p-value"))) +
  ggtitle("viCO coeff., with SVA")+
  annotate("text", x = -0.13, y = -log10(0.12), label = "10% FDR", size = 4, hjust=0)+
  annotate("text", x = -log2(1.5)+0.1, y = 0, label = "<1/1.5 FC", size = 3, hjust=0, vjust=0.5, angle=90)+
  annotate("text", x = -log2(2)+0.1, y = 0, label = "<1/2 FC", size = 3, hjust=0, vjust=0.5,  angle=90)+
  annotate("text", x = -log2(3)+0.1, y = 0, label = "<1/3 FC", size = 3, hjust=0, vjust=0.5, angle=90)+
  annotate("text", x = log2(1.5)+0.1, y = 0, label = ">1.5 FC", size = 3, hjust=0, vjust=0.5, angle=90)+
  annotate("text", x = log2(2)+0.1, y = 0, label = ">2 FC", size = 3, hjust=0, vjust=0.5,  angle=90)+
  annotate("text", x = log2(3)+0.1, y = 0, label = ">3 FC", size = 3, hjust=0, vjust=0.5, angle=90)+
  theme(plot.title = element_text(hjust=0.5, size=16),
        legend.key.size = unit(1.5, "lines"),
        axis.text=element_text(size=14)
  )+
  guides(
    color = guide_legend(override.aes = list(size = 6)),
  )

pdf(file=paste0(dirQCplots,"volcanoPlot_viCO_RG.pdf"), width=8, height=7)
plot(volcanoPlot)
dev.off()


MAplot <- ggplot(data=resFit_conditionvasc, aes(x=AveExpr, y=logFC, col=signif))+
  geom_point(size=1, pch=15, alpha=0.8)+
  theme_bw()+
  scale_color_manual(name="DE genes",values=c("TRUE"="red","FALSE"="grey80"))+
  geom_hline(yintercept=c(-log2(1.5),-log2(2),-log2(3),log2(1.5), log2(2), log2(3)), linetype="dashed", col="grey30")+
  geom_text_repel(data=subset(resFit_conditionvasc, adj.P.Val<0.1), aes(label=gene), col="black")+
  ylab(expression(log[2]*"Fold Change")) +
  xlab(expression(log[2]*"CPM Average expression")) +
  # scale_y_continuous(limits=c(-0.15,0.15), breaks=seq(-0.15,0.15,0.05), labels = scales::label_number())+
  # scale_x_continuous(limits=c(0,12), breaks=seq(0,12,2))+
  ggtitle("vCO coeff., with SVA")+
  annotate("text", y = -log2(1.5)+0.15, x =10.5, label = "<1/1.5x FC", size = 3, hjust=0)+
  annotate("text", y = -log2(2)+0.15, x = 10.5, label = "<1/2x FC", size = 3, hjust=0)+
  annotate("text", y = -log2(3)+0.15, x = 10.5, label = "<1/3x FC", size = 3, hjust=0)+
  annotate("text", y = log2(1.5)+0.15, x = 10.5, label = ">1.5x FC", size = 3, hjust=0)+
  annotate("text", y = log2(2)+0.15, x = 10.5, label = ">2x FC", size = 3, hjust=0)+
  annotate("text", y = log2(3)+0.15, x = 10.5, label = ">3x FC", size = 3, hjust=0)+
  theme(plot.title = element_text(hjust=0.5, size=16),
        legend.key.size = unit(1.5, "lines"),
        axis.text=element_text(size=14)
  )+
  guides(
    color = guide_legend(override.aes = list(size = 6)),
  )


pdf(file=paste0(dirQCplots,"MAPlot_proportionFit_vCO_RG.pdf"), width=8, height=7)
plot(MAplot)
dev.off()


MAplot <- ggplot(data=resFit_conditionvasc_mg, aes(x=AveExpr, y=logFC, col=signif))+
  geom_point(size=1, pch=15, alpha=0.8)+
  theme_bw()+
  scale_color_manual(name="DE genes",values=c("TRUE"="red","FALSE"="grey80"))+
  geom_hline(yintercept=c(-log2(1.5),-log2(2),-log2(3),log2(1.5), log2(2), log2(3)), linetype="dashed", col="grey30")+
  geom_text_repel(data=subset(resFit_conditionvasc_mg, adj.P.Val<0.1), aes(label=gene), col="black")+
  ylab(expression(log[2]*"Fold Change")) +
  xlab(expression(log[2]*"CPM Average expression")) +
  # scale_y_continuous(limits=c(-0.15,0.15), breaks=seq(-0.15,0.15,0.05), labels = scales::label_number())+
  # scale_x_continuous(limits=c(0,12), breaks=seq(0,12,2))+
  ggtitle("viCO coeff., with SVA")+
  annotate("text", y = -log2(1.5)+0.15, x =10.5, label = "<1/1.5x FC", size = 3, hjust=0)+
  annotate("text", y = -log2(2)+0.15, x = 10.5, label = "<1/2x FC", size = 3, hjust=0)+
  annotate("text", y = -log2(3)+0.15, x = 10.5, label = "<1/3x FC", size = 3, hjust=0)+
  annotate("text", y = log2(1.5)+0.15, x = 10.5, label = ">1.5x FC", size = 3, hjust=0)+
  annotate("text", y = log2(2)+0.15, x = 10.5, label = ">2x FC", size = 3, hjust=0)+
  annotate("text", y = log2(3)+0.15, x = 10.5, label = ">3x FC", size = 3, hjust=0)+
  theme(plot.title = element_text(hjust=0.5, size=16),
        legend.key.size = unit(1.5, "lines"),
        axis.text=element_text(size=14)
  )+
  guides(
    color = guide_legend(override.aes = list(size = 6)),
  )


pdf(file=paste0(dirQCplots,"MAPlot_proportionFit_viCO_RG.pdf"),  width=8, height=7)
plot(MAplot)
dev.off()

#######
#######_________________________________________________________________________

### GO enrichment - upregulated and downregulated (for the disease and the two treatment)

geneUniverse <- data.frame(symbol=rownames(se.filt))

library(org.Hs.eg.db)

tabCorr <- AnnotationDbi::select(org.Hs.eg.db, geneUniverse$symbol, "ENTREZID", "SYMBOL")
tabCorr <- tabCorr[!is.na(tabCorr$ENTREZID),]
tabCorr <- tabCorr[!duplicated(tabCorr$SYMBOL),]

geneUniverse$entrezid <- tabCorr[match(geneUniverse$symbol,tabCorr$SYMBOL),]$ENTREZID
geneUniverse <- subset(geneUniverse, !is.na(entrezid))



## function needed

library(GOstats)


GOenrichmentAndReport <- function(geneIds, universeGeneIds,
                                  minSize=3, maxSize=300, minCount=3,
                                  minOddsRatio=1.5, p.value=0.05, highlightGenes=NULL, highlightStr="*%s*",
                                  label="allDE"){
  
  
  GOparams <- new("GOHyperGParams", geneIds=geneIds,
                  universeGeneIds=universeGeneIds,
                  annotation="org.Hs.eg.db",
                  ontology="BP", pvalueCutoff=0.05, conditional=TRUE,
                  minSizeCutoff=3, maxSizeCutoff=300, orCutoff=1.5,
                  testDirection="over")
  
  hgOverGOBP <- hyperGTest(GOparams)
  
  report <- data.frame(GOBPID=as.vector(names(geneIdUniverse(hgOverGOBP))),
                       Pvalue=pvalues(hgOverGOBP),
                       OddsRatio=oddsRatios(hgOverGOBP),
                       ExpCount=expectedCounts(hgOverGOBP),
                       Count=geneCounts(hgOverGOBP),
                       Size=universeCounts(hgOverGOBP),
                       stringsAsFactors=FALSE)
  
  ## discard gene sets that do not meet a minimum and maximum number of genes
  report <- report[report$Size >= minSize & report$Size <= maxSize, , drop=FALSE]
  
  ## discard gene sets that show a p.value>0.05
  report <- report[report$Pvalue < p.value, , drop=FALSE]
  
  ## discard gene sets that do not satisfy the OR cutoff
  report <- report[report$OddsRatio >= minOddsRatio & report$Count >= minCount, , drop=FALSE]
  
  ## apply the maximum-number-of-GO-terms-reported cutoff and sort by odds ratio
  maxReported <- min(nrow(report))
  report <- report[sort(report$OddsRatio, decreasing=TRUE, index.return=TRUE)$ix[1:maxReported], ]
  
  if (dim(report[complete.cases(report),])[1]==0){
    message <- "No GO terms enriched"
    print(message)
    return(message)
  }
  
  ## add the symbol and GO term description, information
  reportGenes <- geneIdsByCategory(hgOverGOBP, report$GOBPID)
  reportGeneSyms <- lapply(reportGenes, annotate::getSYMBOL, "org.Hs.eg.db")
  highlightGeneMask <- lapply(reportGenes, function(x, hgenes, fmt) x %in% hgenes, highlightGenes)
  reportGeneSyms <- mapply(function(genes, mask, fmt) ifelse(mask, sprintf(fmt, genes), genes),
                           reportGeneSyms, highlightGeneMask, MoreArgs=list(fmt=highlightStr), SIMPLIFY=FALSE)
  reportGeneSyms <- sapply(reportGeneSyms, paste, collapse=", ")
  reportGenes <- sapply(reportGenes, paste, collapse=", ")
  report <- cbind(report,
                  Term=sapply(mget(report$GOBPID, GOstats:::GOenv("TERM")), Term),
                  GeneSyms=reportGeneSyms,
                  label=label)
  rownames(report) <- NULL
  
  return(report)
  
}


##

#coefficients_to_runGO <- c("conditionvasc","conditionvasc_mg","donorAlz43cl2","donorPLCG1")
coefficients_to_runGO <- c("conditionvasc","conditionvasc_mg")

cov_goenrich <- sapply(coefficients_to_runGO, function(t) {
  
  print(t)
  
  resFit_coeff <- topTable(fit, coef=t, n=Inf)
  
  upReg_list <- rownames(resFit_coeff[which(resFit_coeff$adj.P.Val<0.1 & resFit_coeff$logFC>0),])
  downReg_list <- rownames(resFit_coeff[which(resFit_coeff$adj.P.Val<0.1 & resFit_coeff$logFC<0),])
  
  #upReg_list <- rownames(resFit_coeff[which(resFit_coeff$logFC>1.5),])
  #downReg_list <- rownames(resFit_coeff[which(resFit_coeff$logFC< -1.5),])
  
  list_up_down <- list(upReg_list,downReg_list)
  names(list_up_down) <- c("Upregulated","Downregulated")
  
  prop_go2 <- sapply(1:length(list_up_down), function(tt){
    
    prop_entrezid <- geneUniverse[match(list_up_down[[tt]], geneUniverse$symbol),]$entrezid
    prop_entrezid <- prop_entrezid[!is.na(prop_entrezid)]
    
    if (length(prop_entrezid)>200){
      
      prop_go <- GOenrichmentAndReport(prop_entrezid, geneUniverse$entrezid,
                                       minSize=5, maxSize=300, minCount=5,
                                       minOddsRatio=1.5, p.value=0.05, highlightGenes=NULL, highlightStr="*%s*",
                                       label=paste0(t,"-", names(list_up_down)[tt]))
    } else {
      
      prop_go <- GOenrichmentAndReport(prop_entrezid, geneUniverse$entrezid,
                                       minSize=3, maxSize=300, minCount=4,
                                       minOddsRatio=1.5, p.value=0.05, highlightGenes=NULL, highlightStr="*%s*",
                                       label=paste0(t,"-", names(list_up_down)[tt]))
    }
    
    if (is.null(dim(prop_go))){
      return(NA)
    } else {
      prop_go$ranking <- 1:dim(prop_go)[1]
      #allres_intersect <- rbind(upreg_go[1:10,], downreg_go[1:10,])
      prop_go$short <- prop_go$Term
      
      return(prop_go)
    }
    
  }, simplify=F)
  
  return(prop_go2)
  
}, simplify=F)

cov_goenrich2 <- unlist(cov_goenrich, recursive = FALSE)
cov_goenrich_high <- cov_goenrich2
cov_goenrich_high <- cov_goenrich_high[sapply(cov_goenrich_high, function(x) !all(is.na(x)))]

cov_goenrich_high <- sapply(cov_goenrich_high, function(x){
  
  if (dim(x)[1]>=10){
    return(x[1:10,])
  } else{
    return(x)
  }
  
}, simplify=F)

cov_goenrich_high <- do.call("rbind", cov_goenrich_high)
rownames(cov_goenrich_high) <- NULL

cov_goenrich_high$short <- cov_goenrich_high$Term
mask_long<- nchar(cov_goenrich_high$Term)>30


cov_goenrich_high$short[mask_long] <- sapply(cov_goenrich_high[mask_long,]$Term, function(x){
  mask <- unlist(sapply(strsplit(x, " "), function(y) cumsum(nchar(y))>30, simplify=F))
  sapply(strsplit(x, " "), function(y) paste0(paste0(y[!mask], collapse=" "),"\n",paste0(y[mask], collapse=" ")))
})

library(tidyverse)
library(tidytext)

cov_goenrich_high <- cov_goenrich_high %>%
  mutate(
    short = reorder_within(short, -ranking, label)
  )

gostats_plot <- ggplot(
  cov_goenrich_high,
  aes(
    x = OddsRatio,
    y = short,
    size = Count,
    col = -log10(Pvalue)
  )
) +
  geom_point() +
  facet_wrap(~label, nrow = 2, scales = "free_y",
             labeller = as_labeller(function(x) {
               x <- gsub("conditionvasc_mg", "</b>viCO</b>", x)
               x <- gsub("conditionvasc", "</b>vCO</b>", x)
               x <- gsub("donorAlz43cl2", "<b>CL2</b>", x)
               x <- gsub("donorPLCG1", "<b>CL3</b>", x)
               x
             })
  ) +
  scale_y_reordered() +
  theme_bw() +
  ylab("") + xlab("Odds ratio") +
  scale_color_viridis_c() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 5),
    axis.title.y = element_text(size = 4.5),
    strip.text = ggtext::element_markdown(size = 7)
  )  +
  ggtitle("GO enrichment per covariate")

pdf(file=paste0(dirQCplots,"GOenrich_DE_RG.pdf"), width=8, height=6)
plot(gostats_plot)
dev.off()




###########################################
########## Heatmaps #######################
##########################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

dirQCplots <- "/scratch/project_2012655/Aurora/plots/RG/"
dir.create(dirQCplots, recursive = TRUE, showWarnings = FALSE)

# 1) Poimi DE-geenit per kontrasti (FDR<0.1), valinnaisesti topN
get_top_up_down_by_effect <- function(tt, fdr = 0.1, n = 10) {
  tt$gene <- rownames(tt)
  
  # 1) pidä vain FDR < 0.1
  tt <- tt %>%
    dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr)
  
  # 2) top UP = suurin positiivinen logFC
  up <- tt %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::arrange(dplyr::desc(logFC), adj.P.Val) %>%
    dplyr::slice_head(n = n)
  
  # 3) top DOWN = negatiivisin logFC
  down <- tt %>%
    dplyr::filter(logFC < 0) %>%
    dplyr::arrange(logFC, adj.P.Val) %>%
    dplyr::slice_head(n = n)
  
  c(up$gene, down$gene)
}

genes_vasc    <- get_top_up_down_by_effect(resFit_conditionvasc,    fdr = 0.1, n = 10)
genes_vasc_mg <- get_top_up_down_by_effect(resFit_conditionvasc_mg, fdr = 0.1, n = 10)

# 2) Tee sample-tason logCPM-matriisista condition-keskiarvot
#    (se.filt sisältää logCPM ja colData(se.filt)$condition)
logCPM <- assays(se.filt)$logCPM  # genes x samples
cond   <- colData(se.filt)$condition
cond   <- factor(cond, levels = c("ctr","vasc","vasc_mg"))

avg_by_condition <- function(mat, cond_factor) {
  stopifnot(ncol(mat) == length(cond_factor))
  out <- sapply(levels(cond_factor), function(cc) {
    rowMeans(mat[, cond_factor == cc, drop = FALSE])
  })
  # out: genes x conditions
  out
}

avg_mat <- avg_by_condition(logCPM, cond)  # genes x (ctr, vasc, vasc_mg)
cond_map <- c("ctr" = "cCO", "vasc" = "vCO", "vasc_mg" = "viCO")
colnames(avg_mat) <- cond_map[colnames(avg_mat)]

# 3) Plot-funktio (genes x conditions)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ggtext)

# ---- 1) Apufunktio: tee heatmap-data yhdelle kontrastille
make_heatmap_df <- function(avg_mat, genes, contrast_label) {
  
  genes <- genes[genes %in% rownames(avg_mat)]
  if (length(genes) == 0) return(NULL)
  
  df <- as.data.frame(avg_mat[genes, c("cCO","vCO","viCO"), drop = FALSE]) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "condition", values_to = "expr")
  
  df$condition <- factor(df$condition, levels = c("cCO","vCO","viCO"))
  
  
  df$contrast <- contrast_label
  df
}


# ---- 2) Rakenna data molemmille kontrasteille
df_vasc <- make_heatmap_df(avg_mat, genes_vasc, "vCO vs cCO")
df_vasc_mg <- make_heatmap_df(avg_mat, genes_vasc_mg, "viCO vs cCO")

heatmap_df <- dplyr::bind_rows(df_vasc, df_vasc_mg) %>%
  dplyr::mutate(
    contrast = factor(contrast, levels = c("vCO vs cCO", "viCO vs cCO")),
    gene_facet = paste(gene, contrast, sep = "___")
  )

levels_vasc <- paste(
  rev(genes_vasc[genes_vasc %in% rownames(avg_mat)]),
  "vCO vs cCO",
  sep = "___"
)

levels_vasc_mg <- paste(
  rev(genes_vasc_mg[genes_vasc_mg %in% rownames(avg_mat)]),
  "viCO vs cCO",
  sep = "___"
)

heatmap_df$gene_facet <- factor(heatmap_df$gene_facet, levels = c(levels_vasc, levels_vasc_mg))

heatmap_df$contrast <- factor(
  heatmap_df$contrast,
  levels = c("vCO vs cCO", "viCO vs cCO")
)




n_down_vasc    <- sum(resFit_conditionvasc[genes_vasc, "logFC"] < 0, na.rm = TRUE)
n_down_vasc_mg <- sum(resFit_conditionvasc_mg[genes_vasc_mg, "logFC"] < 0, na.rm = TRUE)

hline_df <- data.frame(
  contrast = factor(c("vCO vs cCO", "viCO vs cCO"), levels = levels(heatmap_df$contrast)),
  yintercept = c(n_down_vasc + 0.5, n_down_vasc_mg + 0.5)
)


# ---- 3) Piirrä facetoitu heatmap
p_heatmap <- ggplot(heatmap_df, aes(x = condition, y = gene_facet, fill = expr)) +
  geom_tile(color = "black", linewidth = 0.1) +
  geom_hline(
    data = hline_df,
    aes(yintercept = yintercept),
    colour = "white",
    linewidth = 0.75
  ) +
  facet_wrap(~contrast, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  theme_bw() +
  scale_fill_viridis_c() +
  labs(
    title = "DEGs in Radial glia",
    x = "Condition",
    y = "Gene",
    fill = "avg logCPM"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 10)
  )
p_heatmap





library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ggtext)

dirQCplots <- "/scratch/project_2012655/Aurora/plots/RG/"
dir.create(dirQCplots, recursive = TRUE, showWarnings = FALSE)

# 1) Poimi DE-geenit per kontrasti (FDR<0.1), valinnaisesti topN
get_top_up_down_by_effect <- function(tt, fdr = 0.1, n = 10) {
  tt$gene <- rownames(tt)
  
  tt <- tt %>%
    dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr)
  
  up <- tt %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::arrange(dplyr::desc(logFC), adj.P.Val) %>%
    dplyr::slice_head(n = n)
  
  down <- tt %>%
    dplyr::filter(logFC < 0) %>%
    dplyr::arrange(logFC, adj.P.Val) %>%
    dplyr::slice_head(n = n)
  
  # Palauta UP ensin, DOWN perään
  c(up$gene, down$gene)
}

genes_vasc    <- get_top_up_down_by_effect(resFit_conditionvasc,    fdr = 0.1, n = 10)
genes_vasc_mg <- get_top_up_down_by_effect(resFit_conditionvasc_mg, fdr = 0.1, n = 10)

# 2) Condition-keskiarvot logCPM:stä
logCPM <- assays(se.filt)$logCPM  # genes x samples
cond   <- factor(colData(se.filt)$condition, levels = c("ctr","vasc","vasc_mg"))

avg_by_condition <- function(mat, cond_factor) {
  stopifnot(ncol(mat) == length(cond_factor))
  out <- sapply(levels(cond_factor), function(cc) {
    rowMeans(mat[, cond_factor == cc, drop = FALSE])
  })
  out  # genes x conditions
}

avg_mat <- avg_by_condition(logCPM, cond)
cond_map <- c("ctr" = "cCO", "vasc" = "vCO", "vasc_mg" = "viCO")
colnames(avg_mat) <- cond_map[colnames(avg_mat)]

# 3) Tee heatmap-data yhdelle kontrastille (gene merkkijonona)
make_heatmap_df <- function(avg_mat, genes, contrast_label) {
  
  genes_in <- genes[genes %in% rownames(avg_mat)]
  if (length(genes_in) == 0) return(NULL)
  
  df <- as.data.frame(avg_mat[genes_in, c("cCO","vCO","viCO"), drop = FALSE]) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "condition", values_to = "expr")
  
  df$condition <- factor(df$condition, levels = c("cCO","vCO","viCO"))
  df$contrast  <- contrast_label
  df
}

df_vasc    <- make_heatmap_df(avg_mat, genes_vasc,    "vCO vs cCO")
df_vasc_mg <- make_heatmap_df(avg_mat, genes_vasc_mg, "viCO vs cCO")

# Jos jompikumpi on NULL, pudota se pois ennen bind_rows
df_list <- list(df_vasc, df_vasc_mg)
df_list <- df_list[!vapply(df_list, is.null, logical(1))]
heatmap_df <- dplyr::bind_rows(df_list)

# Varmistus: jos tyhjä, ei piirrettävää
stopifnot(nrow(heatmap_df) > 0)

# Facet-kohtainen gene-ID + kontrastin factor
heatmap_df <- heatmap_df %>%
  mutate(
    contrast = factor(contrast, levels = c("vCO vs cCO", "viCO vs cCO")),
    gene_facet = paste(gene, contrast, sep = "___")
  )

# Facet-kohtaiset levelit: rev(UP,DOWN) -> (DOWN,UP) eli UP ylös
levels_vasc <- paste(
  rev(genes_vasc[genes_vasc %in% rownames(avg_mat)]),
  "vCO vs cCO",
  sep = "___"
)

levels_vasc_mg <- paste(
  rev(genes_vasc_mg[genes_vasc_mg %in% rownames(avg_mat)]),
  "viCO vs cCO",
  sep = "___"
)

# Drop NA-levelit (tärkeä jos joku geeni ei päätynyt df:ään)
all_levels <- c(levels_vasc, levels_vasc_mg)
all_levels <- all_levels[!is.na(all_levels)]

heatmap_df$gene_facet <- factor(heatmap_df$gene_facet, levels = all_levels)

# Valkoinen viiva UP/DOWN -rajaan:
# laske down niiden geenien joukosta, jotka oikeasti ovat mukana avg_mat:ssa
genes_vasc_in    <- genes_vasc[genes_vasc %in% rownames(avg_mat)]
genes_vasc_mg_in <- genes_vasc_mg[genes_vasc_mg %in% rownames(avg_mat)]

n_down_vasc <- sum(resFit_conditionvasc[genes_vasc_in, "logFC"] < 0, na.rm = TRUE)
n_down_vasc_mg <- sum(resFit_conditionvasc_mg[genes_vasc_mg_in, "logFC"] < 0, na.rm = TRUE)

hline_df <- data.frame(
  contrast = factor(c("vCO vs cCO", "viCO vs cCO"),
                    levels = levels(heatmap_df$contrast)),
  yintercept = c(n_down_vasc + 0.5, n_down_vasc_mg + 0.5)
)

# 4) Piirrä facetoitu heatmap
p_heatmap <- ggplot(heatmap_df, aes(x = condition, y = gene_facet, fill = expr)) +
  geom_tile(color = "black", linewidth = 0.1) +
  geom_hline(
    data = hline_df,
    aes(yintercept = yintercept),
    colour = "white",
    linewidth = 0.75
  ) +
  facet_wrap(~contrast, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  theme_bw() +
  scale_fill_viridis_c() +
  labs(
    title = "DEGs in Radial glia",
    x = "Condition",
    y = "Gene",
    fill = "avg logCPM"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 10)
  )



pdf(file=paste0(dirQCplots,"heatmap.pdf"), width=8, height=6)
plot(p_heatmap)
dev.off()




###########################################
###   +   Comparison to FindMarkers     ####
###########################################

#PB

pb <- resFit_conditionvasc_mg   # tai resFit_conditionvasc
pb$gene <- rownames(pb)
pb <- pb[!is.na(pb$adj.P.Val), ]

library(Seurat)
library(dplyr)

rg <- subset(obj, subset = sctype_pred == "Radial glia")

# identiteetiksi condition (ctr/vasc/vasc_mg)
Idents(rg) <- rg$condition

fm <- FindMarkers(
  rg,
  ident.1 = "vasc_mg",
  ident.2 = "ctr",
  test.use = "wilcox",    # default
  logfc.threshold = 0,    # älä filtteröi vielä
  min.pct = 0.1,
  return.thresh = 1
)

fm$gene <- rownames(fm)
fm <- fm[!is.na(fm$p_val_adj), ]

if ("avg_log2FC" %in% colnames(fm)) fm$logFC <- fm$avg_log2FC
if ("avg_logFC"  %in% colnames(fm)) fm$logFC <- fm$avg_logFC

#merge

comp <- dplyr::inner_join(
  pb %>% dplyr::select(gene, logFC_pb = logFC, padj_pb = adj.P.Val),
  fm %>% dplyr::select(gene, logFC_fm = logFC, padj_fm = p_val_adj),
  by = "gene"
)

#Top50
topN <- 50

top_pb <- pb %>%
  arrange(adj.P.Val, desc(abs(logFC))) %>%
  slice_head(n = topN) %>% pull(gene)

top_fm <- fm %>%
  arrange(p_val_adj, desc(abs(logFC))) %>%
  slice_head(n = topN) %>% pull(gene)

overlap <- intersect(top_pb, top_fm)
length(overlap)
head(overlap, 20)

# comaprison table 

comp_tbl <- tibble::as_tibble(comp)

report <- comp_tbl %>%
  dplyr::mutate(
    agree_dir = sign(logFC_pb) == sign(logFC_fm),
    rank_pb = rank(padj_pb, ties.method = "min"),
    rank_fm = rank(padj_fm, ties.method = "min")
  ) %>%
  dplyr::arrange(padj_pb) %>%
  dplyr::select(gene, logFC_pb, padj_pb, logFC_fm, padj_fm, agree_dir) %>%
  dplyr::slice_head(n = 50)

utils::write.csv(
  report,
  file = file.path(dirQCplots, "compare_pseudobulk_vs_FindMarkers_RG_vascmg_vs_ctr.csv"),
  row.names = FALSE
)



########################
## 2. GSVA 3.script ####
########################

# Huom pseudobulking kannattaa tehdä eri lailla GSVA:han kuin tähän -> 3. scriptin 0 ja 3 vaiheet









##########################
## 3. FindMaerkers to SMC ##
#########################

smc <- subset(obj, subset = sctype_pred == "Smooth muscle")
stopifnot(all(c("vasc", "vasc_mg") %in% unique(as.character(smc$condition))))
Idents(smc) <- smc$condition

fm <- FindMarkers(
  smc,
  ident.1 = "vasc",
  ident.2 = "vasc_mg",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1,
  return.thresh = 1
)

fm2 <- fm

# gene-nimi
fm2$gene <- rownames(fm2)

# tee logFC-sarake (Seurat version mukaan)
if ("avg_log2FC" %in% colnames(fm2)) fm2$logFC <- fm2$avg_log2FC
if ("avg_logFC"  %in% colnames(fm2)) fm2$logFC <- fm2$avg_logFC

# tee adj.P.Val-sarake (limma-tyyliin)
fm2$adj.P.Val <- fm2$p_val_adj

# signif (limma-tyyliin)
fm2$signif <- ifelse(!is.na(fm2$adj.P.Val) & fm2$adj.P.Val < 0.1, TRUE, FALSE)


volcanoPlot <- ggplot(data=fm2, aes(x=logFC, y=-log10(adj.P.Val), col=signif))+
  geom_point(size=1, pch=15)+
  theme_bw()+
  scale_color_manual(name="DE genes",values=c("TRUE"="red","FALSE"="grey80"))+
  geom_hline(yintercept=-log10(0.1), linetype="dashed", col="grey30")+
  geom_vline(xintercept=c(-log2(1.5),-log2(2),-log2(3), log2(1.5), log2(2), log2(3), log2(4)),
             linetype="dashed", col="grey30")+
  geom_text_repel(
    data=subset(fm2, adj.P.Val < 0.1),
    aes(label=gene),
    col="black",
    size=3,
    max.overlaps = 30
  )+
  xlab(expression(log[2]*"Fold Change")) +
  ylab(expression(-log[10]("adj. p-value"))) +
  ggtitle("SMC: vco vs vico (Wilcox)")+
  annotate("text", x = -0.13, y = -log10(0.12), label = "10% FDR", size = 4, hjust=0)+
  annotate("text", x = -log2(1.5)+0.1, y = 0, label = "<1/1.5 FC", size = 3, hjust=0, vjust=0.5, angle=90)+
  annotate("text", x = -log2(2)+0.1,   y = 0, label = "<1/2 FC",   size = 3, hjust=0, vjust=0.5, angle=90)+
  annotate("text", x = -log2(3)+0.1,   y = 0, label = "<1/3 FC",   size = 3, hjust=0, vjust=0.5, angle=90)+
  annotate("text", x =  log2(1.5)+0.1, y = 0, label = ">1.5 FC",   size = 3, hjust=0, vjust=0.5, angle=90)+
  annotate("text", x =  log2(2)+0.1,   y = 0, label = ">2 FC",     size = 3, hjust=0, vjust=0.5, angle=90)+
  annotate("text", x =  log2(3)+0.1,   y = 0, label = ">3 FC",     size = 3, hjust=0, vjust=0.5, angle=90)+
  theme(plot.title = element_text(hjust=0.5, size=16),
        legend.key.size = unit(1.5, "lines"),
        axis.text=element_text(size=14)
  )+
  guides(color = guide_legend(override.aes = list(size = 6)))

dirQCplots <- "/scratch/project_2012655/Aurora/plots/SMC_FB_EC/"

pdf(file=paste0(dirQCplots, "volcano_SMC_vasc_vs_vasc_mg_wilcox.pdf"), width=8, height=7)
plot(volcanoPlot)
dev.off()


############################################
## 3. GOstats GO BP enrichment from FindMarkers
## (Seurat Wilcox: vasc vs vasc_mg)
############################################

library(org.Hs.eg.db)
library(AnnotationDbi)
library(GOstats)
library(tidyverse)
library(tidytext)
library(ggtext)
library(ggplot2)

## 0) Standardisoi FindMarkers-output "limma-tyyliin"
fm2 <- fm
fm2$gene <- rownames(fm2)

# logFC sarake (Seurat-versiosta riippuen)
if ("avg_log2FC" %in% colnames(fm2)) fm2$logFC <- fm2$avg_log2FC
if ("avg_logFC"  %in% colnames(fm2)) fm2$logFC <- fm2$avg_logFC
stopifnot("logFC" %in% colnames(fm2))

# adj.P.Val sarake
stopifnot("p_val_adj" %in% colnames(fm2))
fm2$adj.P.Val <- fm2$p_val_adj

## 1) geneUniverse (taustageenit) = kaikki FindMarkersissa testatut geenit
geneUniverse <- data.frame(symbol = fm2$gene)

tabCorr <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = geneUniverse$symbol,
  columns  = c("ENTREZID", "SYMBOL"),
  keytype  = "SYMBOL"
)
tabCorr <- tabCorr[!is.na(tabCorr$ENTREZID), ]
tabCorr <- tabCorr[!duplicated(tabCorr$SYMBOL), ]

geneUniverse$entrezid <- tabCorr[match(geneUniverse$symbol, tabCorr$SYMBOL), "ENTREZID"]
geneUniverse <- subset(geneUniverse, !is.na(entrezid))

## 2) GOstats-funktio (sun koodi)
GOenrichmentAndReport <- function(geneIds, universeGeneIds,
                                  minSize=3, maxSize=300, minCount=3,
                                  minOddsRatio=1.5, p.value=0.05, highlightGenes=NULL, highlightStr="*%s*",
                                  label="allDE"){
  
  GOparams <- new("GOHyperGParams", geneIds=geneIds,
                  universeGeneIds=universeGeneIds,
                  annotation="org.Hs.eg.db",
                  ontology="BP", pvalueCutoff=p.value, conditional=TRUE,
                  minSizeCutoff=minSize, maxSizeCutoff=maxSize,
                  testDirection="over")
  
  hgOverGOBP <- hyperGTest(GOparams)
  
  report <- data.frame(GOBPID=as.vector(names(geneIdUniverse(hgOverGOBP))),
                       Pvalue=pvalues(hgOverGOBP),
                       OddsRatio=oddsRatios(hgOverGOBP),
                       ExpCount=expectedCounts(hgOverGOBP),
                       Count=geneCounts(hgOverGOBP),
                       Size=universeCounts(hgOverGOBP),
                       stringsAsFactors=FALSE)
  
  report <- report[report$Size >= minSize & report$Size <= maxSize, , drop=FALSE]
  report <- report[report$Pvalue < p.value, , drop=FALSE]
  report <- report[report$OddsRatio >= minOddsRatio & report$Count >= minCount, , drop=FALSE]
  
  if (nrow(report[complete.cases(report), , drop=FALSE]) == 0){
    message("No GO terms enriched: ", label)
    return(NULL)
  }
  
  reportGenes <- geneIdsByCategory(hgOverGOBP, report$GOBPID)
  reportGeneSyms <- lapply(reportGenes, annotate::getSYMBOL, "org.Hs.eg.db")
  
  highlightGeneMask <- lapply(reportGenes, function(x, hgenes) x %in% hgenes, highlightGenes)
  reportGeneSyms <- mapply(function(genes, mask, fmt) ifelse(mask, sprintf(fmt, genes), genes),
                           reportGeneSyms, highlightGeneMask,
                           MoreArgs=list(fmt=highlightStr), SIMPLIFY=FALSE)
  
  reportGeneSyms <- sapply(reportGeneSyms, paste, collapse=", ")
  report <- cbind(report,
                  Term = sapply(mget(report$GOBPID, GOstats:::GOenv("TERM")), Term),
                  GeneSyms = reportGeneSyms,
                  label = label)
  rownames(report) <- NULL
  report
}

## 3) Up/Down geenilistat FindMarkersista
fdr_cut <- 0.1  # sama kuin sun muissa
upReg_list <- fm2$gene[!is.na(fm2$adj.P.Val) & fm2$adj.P.Val < fdr_cut & fm2$logFC > 0]
downReg_list <- fm2$gene[!is.na(fm2$adj.P.Val) & fm2$adj.P.Val < fdr_cut & fm2$logFC < 0]

list_up_down <- list(Upregulated = upReg_list, Downregulated = downReg_list)

## 4) Aja GO Up ja Down (samalla logiikalla kuin sulla)
go_res <- lapply(names(list_up_down), function(tt){
  
  prop_entrezid <- geneUniverse$entrezid[match(list_up_down[[tt]], geneUniverse$symbol)]
  prop_entrezid <- prop_entrezid[!is.na(prop_entrezid)]
  
  if (length(prop_entrezid) == 0) return(NULL)
  
  if (length(prop_entrezid) > 200){
    prop_go <- GOenrichmentAndReport(prop_entrezid, geneUniverse$entrezid,
                                     minSize=5, maxSize=300, minCount=5,
                                     minOddsRatio=1.5, p.value=0.05,
                                     highlightGenes=NULL, highlightStr="*%s*",
                                     label=paste0("vasc_vs_vasc_mg-", tt))
  } else {
    prop_go <- GOenrichmentAndReport(prop_entrezid, geneUniverse$entrezid,
                                     minSize=3, maxSize=300, minCount=4,
                                     minOddsRatio=1.5, p.value=0.05,
                                     highlightGenes=NULL, highlightStr="*%s*",
                                     label=paste0("vasc_vs_vasc_mg-", tt))
  }
  
  if (is.null(prop_go)) return(NULL)
  
  prop_go$ranking <- seq_len(nrow(prop_go))
  prop_go$short <- prop_go$Term
  prop_go
})
names(go_res) <- names(list_up_down)

go_res <- go_res[!sapply(go_res, is.null)]



## 6) Tee dotplot top10 / suunta (sun plot-tyyli)
if (length(go_res) == 0) {
  message("No GO results to plot.")
} else {
  
  cov_goenrich_high <- lapply(go_res, function(x){
    if (nrow(x) >= 10) x[1:10, , drop=FALSE] else x
  })
  cov_goenrich_high <- do.call("rbind", cov_goenrich_high)
  rownames(cov_goenrich_high) <- NULL
  
  cov_goenrich_high$short <- cov_goenrich_high$Term
  mask_long <- nchar(cov_goenrich_high$Term) > 30
  
  cov_goenrich_high$short[mask_long] <- sapply(cov_goenrich_high[mask_long, ]$Term, function(x){
    mask <- unlist(sapply(strsplit(x, " "), function(y) cumsum(nchar(y)) > 30, simplify=FALSE))
    sapply(strsplit(x, " "), function(y) paste0(paste0(y[!mask], collapse=" "), "\n", paste0(y[mask], collapse=" ")))
  })
  
  cov_goenrich_high <- cov_goenrich_high %>%
    mutate(short = reorder_within(short, -ranking, label))
  
  gostats_plot <- ggplot(
    cov_goenrich_high,
    aes(x = OddsRatio, y = short, size = Count, col = -log10(Pvalue))
  ) +
    geom_point() +
    facet_wrap(~label, nrow = 2, scales = "free_y") +
    scale_y_reordered() +
    theme_bw() +
    ylab("") + xlab("Odds ratio") +
    scale_color_viridis_c() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 5),
      axis.title.y = element_text(size = 4.5),
      strip.text = ggtext::element_markdown(size = 7)
    ) +
    ggtitle("GO enrichment (FindMarkers DE genes)")
  
  pdf(file=paste0(dirQCplots,"GOenrich_DE_fromFindMarkers.pdf"), width=8, height=6)
  plot(gostats_plot)
  dev.off()
}



#### 

############################################################
## 4. RG FindMarkers per donor (vasc_mg vs ctr) + intersect
## + donor-donor logFC scatter plots (1-2, 1-3, 2-3)
############################################################

# --- libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# --- user settings
dirQCplots <- "/scratch/project_2012655/Aurora/plots/RG/"
fdr_cut    <- 0.1
min_cells  <- 20

# --- 0) subset RG + set Idents
rg <- subset(obj, subset = sctype_pred == "Radial glia")
Idents(rg) <- rg$condition

# donors (your 3 donors)
donors <- c("Alz37cl2", "Alz43cl2", "PLCG1")
donors <- donors[donors %in% unique(rg$donor)]
stopifnot(length(donors) >= 3)

message("Cells per donor x condition (RG):")
print(table(rg$donor, rg$condition))

# --- 1.helper functions: FindMarkers for one donor
run_fm_one_donor <- function(seurat_obj, donor_id,
                             ident.1 = "vasc_mg", ident.2 = "ctr",
                             min_cells = 20,
                             test.use = "wilcox",
                             logfc.threshold = 0,
                             min.pct = 0.1,
                             return.thresh = 1) {
  
  x <- subset(seurat_obj, subset = donor == donor_id)
  
  tab <- table(Idents(x))
  if (!all(c(ident.1, ident.2) %in% names(tab))) {
    message("SKIP ", donor_id, " (missing groups: ", ident.1, " or ", ident.2, ")")
    return(NULL)
  }
  if (tab[[ident.1]] < min_cells || tab[[ident.2]] < min_cells) {
    message("SKIP ", donor_id, " (too few cells: ",
            ident.1, "=", tab[[ident.1]], ", ", ident.2, "=", tab[[ident.2]], ")")
    return(NULL)
  }
  
  fm <- FindMarkers(
    x,
    ident.1 = ident.1,
    ident.2 = ident.2,
    test.use = test.use,
    logfc.threshold = logfc.threshold,
    min.pct = min.pct,
    return.thresh = return.thresh
  )
  
  fm$gene <- rownames(fm)
  
  # unify logFC column name
  if ("avg_log2FC" %in% colnames(fm)) fm$logFC <- fm$avg_log2FC
  if ("avg_logFC"  %in% colnames(fm)) fm$logFC <- fm$avg_logFC
  
  fm$donor <- donor_id
  fm
}

# --- 2. helper: significant genes (FDR)
get_sig_genes <- function(fm, fdr = 0.1, lfc_min = 0) {
  fm %>%
    dplyr::filter(!is.na(p_val_adj)) %>%
    dplyr::filter(p_val_adj < fdr) %>%
    dplyr::filter(abs(logFC) > lfc_min) %>%
    dplyr::pull(gene) %>%
    unique()
}

# --- 3.  helper: build wide logFC table for a gene set, for selected donors
build_lfc_wide <- function(fm_list, donors_use, genes_use) {
  
  long <- dplyr::bind_rows(lapply(donors_use, function(d) {
    fm_list[[d]] %>%
      dplyr::select(gene, logFC) %>%
      dplyr::filter(gene %in% genes_use) %>%
      dplyr::mutate(donor = d)
  })) %>%
    dplyr::distinct(gene, donor, .keep_all = TRUE)
  
  wide <- long %>%
    tidyr::pivot_wider(names_from = donor, values_from = logFC, values_fn = mean)
  
  # ensure numeric
  dcols <- setdiff(colnames(wide), "gene")
  wide[dcols] <- lapply(wide[dcols], as.numeric)
  
  wide
}

# --- 4. helper: plot scatter
plot_scatter <- function(df_wide, d1, d2, title, subtitle_prefix = "vasc_mg vs ctr") {
  
  tmp <- df_wide %>%
    dplyr::filter(!is.na(.data[[d1]]), !is.na(.data[[d2]]))
  
  if (nrow(tmp) < 3) {
    return(ggplot() + theme_void() + ggtitle(title) +
             labs(subtitle = paste0("Too few genes (n=", nrow(tmp), ")")))
  }
  
  cor_val <- cor(tmp[[d1]], tmp[[d2]], use = "pairwise.complete.obs")
  
  ggplot(tmp, aes(x = .data[[d1]], y = .data[[d2]])) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    theme_bw() +
    labs(
      title = title,
      subtitle = paste0(subtitle_prefix,
                        " | genes = ", nrow(tmp),
                        " | Pearson r = ", round(cor_val, 2)),
      x = paste0("logFC (", d1, ")"),
      y = paste0("logFC (", d2, ")")
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

# --- 1) run FindMarkers for each donor
fm_list <- lapply(donors, function(d) run_fm_one_donor(
  rg, d, ident.1 = "vasc_mg", ident.2 = "ctr", min_cells = min_cells
))
names(fm_list) <- donors
fm_list <- fm_list[!sapply(fm_list, is.null)]

# Saving
for (d in names(fm_list)) {
  
  fm <- fm_list[[d]]
  
  # CSV (helppo lukea, jakaa)
  write.csv(
    fm,
    file = file.path(
      dirQCplots,
      paste0("RG_FindMarkers_", d, "_vasc_mg_vs_ctr_Wilcoxon.csv")
    ),
    row.names = FALSE
  )
  
  # RDS (täydellinen, nopea ladata R:ään)
  saveRDS(
    fm,
    file = file.path(
      dirQCplots,
      paste0("RG_FindMarkers_", d, "_vasc_mg_vs_ctr_Wilcoxon.rds")
    )
  )
}
stopifnot(length(fm_list) >= 3)



# --- 2) significant genes per donor
sig_list <- lapply(fm_list, get_sig_genes, fdr = fdr_cut, lfc_min = 0)
n_de <- sapply(sig_list, length)
message("N DE genes per donor (FDR < ", fdr_cut, "):")
print(n_de)

# --- 3) make 3 pairwise scatters (genes = intersect of the pair)
pairs <- combn(names(fm_list), 2, simplify = FALSE)

pair_plots <- lapply(pairs, function(p) {
  d1 <- p[1]; d2 <- p[2]
  genes_pair <- union(sig_list[[d1]], sig_list[[d2]])
  
  df_wide <- build_lfc_wide(fm_list, donors_use = c(d1, d2), genes_use = genes_pair)
  
  plot_scatter(
    df_wide, d1, d2,
    title = paste0("Pairwise intersect: ", d1, " vs ", d2)
  )
})
names(pair_plots) <- sapply(pairs, function(p) paste(p, collapse = "_vs_"))

# --- 4) all-3 intersect scatter (genes = intersect across all 3 donors)
d_all <- names(fm_list)[1:3]  # three donors
genes_all3 <- Reduce(intersect, sig_list[d_all])

# choose axes for the "all donors intersect" plot (use donor1 vs donor2)
d1_all <- d_all[1]
d2_all <- d_all[2]

df_all3 <- build_lfc_wide(fm_list, donors_use = c(d1_all, d2_all), genes_use = genes_all3)

p_all3 <- plot_scatter(
  df_all3, d1_all, d2_all,
  title = paste0("ALL-3 intersect genes (", paste(d_all, collapse = ", "), ")"),
  subtitle_prefix = "vasc_mg vs ctr | genes DE in ALL 3 donors"
)

# --- 5) combine: 3 pairwise + 1 all3 = 4 plots
plots_4 <- c(pair_plots, list(ALL3 = p_all3))

p_out <- wrap_plots(plots_4, ncol = 2) +
  plot_annotation(
    title = paste0("RG FindMarkers logFC consistency (Wilcoxon)\n",
                   "Pairwise intersect per donor-pair + all-3 intersect (FDR < ", fdr_cut, ")"),
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
  )

# save
pdf(
  file.path(dirQCplots, paste0("RG_FindMarkers_logFC_4scatters_pairwisePlusAll3_vasc_mg_vs_ctr_FDR", fdr_cut, ".pdf")),
  width = 12,
  height = 9
)
print(p_out)
dev.off()

p_out






# versio 2
# --- libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

# --- user settings
dirQCplots <- "/scratch/project_2012655/Aurora/plots/RG/"
fdr_cut    <- 0.1
min_cells  <- 20

# --- 0) subset RG + set Idents
rg <- subset(obj, subset = sctype_pred == "Radial glia")
Idents(rg) <- rg$condition

# donors (your 3 donors)
donors <- c("Alz37cl2", "Alz43cl2", "PLCG1")
donors <- donors[donors %in% unique(rg$donor)]
stopifnot(length(donors) >= 3)

message("Cells per donor x condition (RG):")
print(table(rg$donor, rg$condition))

# --- 1) helper: FindMarkers for one donor /vico vs cco
run_fm_one_donor <- function(seurat_obj, donor_id,
                             ident.1 = "vasc_mg", ident.2 = "ctr",
                             min_cells = 20,
                             test.use = "wilcox",
                             logfc.threshold = 0,
                             min.pct = 0.1,
                             return.thresh = 1) {
  
  x <- subset(seurat_obj, subset = donor == donor_id)
  
  tab <- table(Idents(x))
  if (!all(c(ident.1, ident.2) %in% names(tab))) {
    message("SKIP ", donor_id, " (missing groups: ", ident.1, " or ", ident.2, ")")
    return(NULL)
  }
  if (tab[[ident.1]] < min_cells || tab[[ident.2]] < min_cells) {
    message("SKIP ", donor_id, " (too few cells: ",
            ident.1, "=", tab[[ident.1]], ", ", ident.2, "=", tab[[ident.2]], ")")
    return(NULL)
  }
  
  fm <- FindMarkers(
    x,
    ident.1 = ident.1,
    ident.2 = ident.2,
    test.use = test.use,
    logfc.threshold = logfc.threshold,
    min.pct = min.pct,
    return.thresh = return.thresh
  )
  
  fm$gene <- rownames(fm)
  
  # unify logFC column name
  if ("avg_log2FC" %in% colnames(fm)) fm$logFC <- fm$avg_log2FC
  if ("avg_logFC"  %in% colnames(fm)) fm$logFC <- fm$avg_logFC
  
  fm$donor <- donor_id
  fm
}

# --- 2) helper: significant genes (FDR)
get_sig_genes <- function(fm, fdr = 0.1, lfc_min = 0) {
  fm %>%
    dplyr::filter(!is.na(p_val_adj)) %>%
    dplyr::filter(p_val_adj < fdr) %>%
    dplyr::filter(abs(logFC) > lfc_min) %>%
    dplyr::pull(gene) %>%
    unique()
}

# --- 3) run FindMarkers for each donor
fm_list <- lapply(donors, function(d) run_fm_one_donor(
  rg, d, ident.1 = "vasc_mg", ident.2 = "ctr", min_cells = min_cells
))
names(fm_list) <- donors
fm_list <- fm_list[!sapply(fm_list, is.null)]
stopifnot(length(fm_list) >= 3)


# --- 5) significant genes per donor
sig_list <- lapply(fm_list, get_sig_genes, fdr = fdr_cut, lfc_min = 0)
n_de <- sapply(sig_list, length)
message("N DE genes per donor (FDR < ", fdr_cut, "):")
print(n_de)

# --- 6) Build one dataframe for ALL donor-pairs (UNION + 3-color groups)
donor_labels <- c(
  "Alz37cl2" = "CL1",
  "Alz43cl2" = "CL2",
  "PLCG1"    = "CL3"
)

# sanity: ensure mapping exists for all donors used
stopifnot(all(names(fm_list) %in% names(donor_labels)))

build_pairwise_df <- function(fm_list, sig_list, donor_labels) {
  
  donors <- as.character(names(fm_list))
  stopifnot(length(donors) >= 2)
  
  pair_mat <- t(combn(donors, 2))  # rows: (d1,d2)
  
  df_all <- dplyr::bind_rows(apply(pair_mat, 1, function(row) {
    
    d1 <- row[[1]]
    d2 <- row[[2]]
    
    genes_union  <- union(sig_list[[d1]], sig_list[[d2]])
    genes_shared <- intersect(sig_list[[d1]], sig_list[[d2]])
    genes_d1only <- setdiff(sig_list[[d1]], sig_list[[d2]])
    genes_d2only <- setdiff(sig_list[[d2]], sig_list[[d1]])
    
    l1 <- fm_list[[d1]] %>% dplyr::select(gene, logFC) %>% dplyr::distinct(gene, .keep_all = TRUE)
    l2 <- fm_list[[d2]] %>% dplyr::select(gene, logFC) %>% dplyr::distinct(gene, .keep_all = TRUE)
    
    dplyr::full_join(l1, l2, by = "gene", suffix = c("_d1", "_d2")) %>%
      dplyr::filter(gene %in% genes_union) %>%
      dplyr::rename(logFC_d1 = logFC_d1, logFC_d2 = logFC_d2) %>%
      dplyr::mutate(
        pair = paste0(donor_labels[[d1]], " vs ", donor_labels[[d2]]),
        group = dplyr::case_when(
          gene %in% genes_shared ~ "shared",
          gene %in% genes_d1only ~ paste0(donor_labels[[d1]], "_only"),
          gene %in% genes_d2only ~ paste0(donor_labels[[d2]], "_only"),
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(!is.na(logFC_d1), !is.na(logFC_d2), !is.na(group))
    
  }))
  
  df_all
}

# ✅ IMPORTANT: pass donor_labels as 3rd argument
df_pairs <- build_pairwise_df(fm_list, sig_list, donor_labels)



# --- 8) Plot facets: donor-pair shown in facet strip title
p_facets <- ggplot(df_pairs, aes(x = logFC_d1, y = logFC_d2, color = group)) +
  geom_point(alpha = 0.8, size = 1.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  facet_wrap(~ pair, scales = "free") +
  theme_bw() +
  scale_color_brewer(palette = "Set1", name = "DE group") +
  labs(
    title = paste0("RG donor-pair logFC consistency (viCO vs cCO)"),
    x = "logFC (donor 1 in pair)",
    y = "logFC (donor 2 in pair)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )
p_facets

# --- 9) Save
out_pdf <- file.path(dirQCplots, paste0("RG_FindMarkers_logFC_facets_pairwiseUNION_3colors_FDR", fdr_cut, ".pdf"))
ggsave(filename = out_pdf, plot = p_facets, width = 12, height = 5)



# Heatmap per donor --> LOG FOLD changes GENES --> per donor --> expressio in cco (donor)














# 5. HEATMAP genes RG

#==============================

suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(viridis)
  library(ggh4x)
  library(scales)
})

# ---- paths
dirQCplots <- "/scratch/project_2012655/Aurora/plots/RG/"
dir.create(dirQCplots, recursive = TRUE, showWarnings = FALSE)
out_pdf <- file.path(dirQCplots, "heatmap_RG_markers_logCPM_facetsCL_xConditionOnly.pdf")

# ---- genes
genes_use <- c(
  "SOX9","NFIA","NFIB","HES1","HES5",
  "ASCL1","NEUROG1","NEUROG2","DLL1","HES6",
  "DCX","PAX6"
)

# ---- condition ordering + display labels
cond_levels <- c("ctr","vasc","vasc_mg")
cond_disp   <- c("ctr"="cCO", "vasc"="vCO", "vasc_mg"="viCO")

# ===============================
# 0) checks
# ===============================
stopifnot(exists("pseudo_obj_sub_RG"), inherits(pseudo_obj_sub_RG, "Seurat"))

# ===============================
# 1) Pseudobulk counts + metadata (RG)
# ===============================
meta <- pseudo_obj_sub_RG@meta.data
req_cols <- c("condition","orig.ident")
miss <- setdiff(req_cols, colnames(meta))
if (length(miss) > 0) stop("Missing meta columns in pseudo_obj_sub_RG@meta.data: ", paste(miss, collapse = ", "))

counts <- as.matrix(GetAssayData(pseudo_obj_sub_RG, assay = "RNA", slot = "counts"))
if (!all(colnames(counts) %in% rownames(meta))) stop("counts colnames do not match rownames(meta).")
meta <- meta[colnames(counts), , drop = FALSE]

# condition factor
meta$condition <- factor(as.character(meta$condition), levels = cond_levels)

# sample_id = part before "_" (AggregateExpression often makes "SAMPLE_Celltype")
meta$sample_id <- sub("_(.*)$", "", as.character(meta$orig.ident))

# cell line label (display only)
meta$cellline <- dplyr::case_when(
  grepl("Alz37cl2", meta$sample_id) ~ "CL1",
  grepl("Alz43cl2", meta$sample_id) ~ "CL2",
  grepl("PLCG1",    meta$sample_id) ~ "CL3",
  TRUE ~ meta$sample_id
)
meta$cellline <- factor(meta$cellline, levels = c("CL1","CL2","CL3"))

# ===============================
# 2) EXPRESSION = logCPM from counts (edgeR)
# ===============================
y <- DGEList(counts = counts)
y <- calcNormFactors(y)
logCPM <- cpm(y, log = TRUE, prior.count = 1)  # genes x pseudobulk samples

# ===============================
# 3) Long df: gene x (cellline, condition)
# ===============================
genes_present <- intersect(genes_use, rownames(logCPM))
if (length(genes_present) == 0) stop("None of the requested genes found in logCPM rownames.")

heat_df <- logCPM[genes_present, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "pb_sample", values_to = "expr") %>%
  mutate(
    condition_raw = as.character(meta$condition[match(pb_sample, rownames(meta))]),
    condition     = factor(cond_disp[condition_raw], levels = c("cCO","vCO","viCO")),
    cellline      = meta$cellline[match(pb_sample, rownames(meta))],
    gene          = factor(gene, levels = rev(genes_use))
  )


# Clip color scale (recommended; does not change values)
lims <- quantile(heat_df$expr, probs = c(0.05, 0.95), na.rm = TRUE)


# ===============================
# ADD: RG cell counts per (cellline x condition)
# ===============================
rg_counts <- obj_sub@meta.data %>%
  filter(sctype_pred == "Radial glia") %>%
  mutate(
    sample_id = as.character(orig.ident),
    cellline = case_when(
      grepl("Alz37cl2", sample_id) ~ "CL1",
      grepl("Alz43cl2", sample_id) ~ "CL2",
      grepl("PLCG1",    sample_id) ~ "CL3",
      TRUE ~ NA_character_
    ),
    condition = factor(condition, levels = cond_levels),
    condition = cond_disp[as.character(condition)]
  ) %>%
  group_by(cellline, condition) %>%
  summarise(n_RG = n(), .groups = "drop")

# ===============================
# 4) Plot: TOP = CL facets, BOTTOM = condition only
# ===============================
p <- ggplot(heat_df, aes(x = condition, y = gene, fill = expr)) +
  geom_tile(color = "black", linewidth = 0.1) +
  ggh4x::facet_grid2(
    . ~ cellline,
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_fill_viridis_c(limits = lims, oob = scales::squish) +
  theme_bw() +
  labs(
    title = "RG marker expression",
    x = "",
    y = "",
    fill = "logCPM"
  ) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    strip.text  = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    panel.grid  = element_blank()
  )
p <- p +
  geom_text(
    data = rg_counts,
    aes(x = condition, y = -0.2, label = paste0("n=", n_RG)),
    inherit.aes = FALSE,
    size = 3,
    vjust = 1
  )
print(p)
p <- p +
  geom_text(
    data = rg_counts,
    aes(x = condition, y = 0.2, label = paste0("n=", n_RG)),
    inherit.aes = FALSE,
    size = 3,
    vjust = 1
  )
ggsave(out_pdf, p, device = cairo_pdf, width = 7.0, height = 4.8)
print(p)
message("Saved: ", out_pdf)

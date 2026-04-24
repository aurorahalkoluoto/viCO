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

#################
#Fig 3E-F, S4A##
#################

obj2 <- readRDS()


ctypes_to_discard <- table(obj$orig.ident, obj$sctype_pred)>=10
ctypes_to_discard <- names(which(colSums(ctypes_to_discard)<dim(ctypes_to_discard)[1]))

obj_sub <- obj[,!obj$sctype_pred %in% ctypes_to_discard]

stopifnot(all(table(obj_sub$orig.ident, obj_sub$sctype_pred)>=10))

## Pseudobulk by sample

pseudo_obj_sub <- AggregateExpression(
  obj_sub,
  group.by = c("orig.ident","sctype_pred"),  
  assays = "RNA",
  slot = "counts",  
  return.seurat=T
)


pseudo_obj_sub <- FindVariableFeatures(pseudo_obj_sub)




pseudo_obj_sub$old.ident <- sapply(strsplit(pseudo_obj_sub$orig.ident,"_"), function(x) x[1])
obj$orig.ident <- gsub("_","-", obj$orig.ident)

pseudo_obj_sub$batch <- unname(obj$batch[match(pseudo_obj_sub$old.ident, obj$orig.ident)])
pseudo_obj_sub$condition <- unname(obj$condition[match(pseudo_obj_sub$old.ident, obj$orig.ident)])
pseudo_obj_sub$donor <- unname(obj$donor[match(pseudo_obj_sub$old.ident, obj$orig.ident)])




#  Radial glia (DEGs+GO)


pseudo_obj_sub_RG <- subset(pseudo_obj_sub, sctype_pred=="Radial glia")


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


dge.filt.noTMM <- dge.filt

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


df_summary <- as.data.frame(summary(res))
colnames(df_summary) <- c("Direction","Covariate","numDE")

### Plot per covariate

resFit_conditionvasc <- topTable(fit, coef="conditionvasc", n=Inf)
resFit_conditionvasc_mg <- topTable(fit, coef="conditionvasc_mg", n=Inf)


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



coefficients_to_runGO <- c("conditionvasc","conditionvasc_mg")

cov_goenrich <- sapply(coefficients_to_runGO, function(t) {
  
  print(t)
  
  resFit_coeff <- topTable(fit, coef=t, n=Inf)
  
  upReg_list <- rownames(resFit_coeff[which(resFit_coeff$adj.P.Val<0.1 & resFit_coeff$logFC>0),])
  downReg_list <- rownames(resFit_coeff[which(resFit_coeff$adj.P.Val<0.1 & resFit_coeff$logFC<0),])
  
  
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

###############
## FIGURE 3G ##
###############


# HEATMAP genes RG

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

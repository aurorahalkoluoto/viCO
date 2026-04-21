

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
dirQCplots <- "/scratch/project_2012655/Aurora/plots/all/"


###########################
## Pseudobulk per sample ## 
###########################

## Sample: Each donor in the three conditions (9 samples). Cell types are ignored
## The sample is treaten like a juice with all the fruits/cell types

## Check that at least all samples have 10 cells
stopifnot(!any(names(which(table(obj$orig.ident)<10))))

pseudo_obj <- AggregateExpression(
  obj,
  group.by = "orig.ident",  # Replace with your grouping variable (e.g., sample, condition)
  assays = "RNA",
  slot = "counts",  # Use log-normalized expression instead of raw counts
  return.seurat=T
)


## Include information about batch, condition and donor

obj$orig.ident <- gsub("_","-",obj$orig.ident)

pseudo_obj$batch <- unname(obj$batch[match(pseudo_obj$orig.ident, obj$orig.ident)])
pseudo_obj$condition <- unname(obj$condition[match(pseudo_obj$orig.ident, obj$orig.ident)])
pseudo_obj$donor <- unname(obj$donor[match(pseudo_obj$orig.ident, obj$orig.ident)])

## Find highly variable genes
pseudo_obj <- FindVariableFeatures(pseudo_obj)


library(SummarizedExperiment)
library(edgeR)
library(limma)
library(sva)

#### differential expression

se <- SummarizedExperiment(assays=list(counts=as.matrix(pseudo_obj@assays$RNA["counts"])),
                           colData=pseudo_obj@meta.data)

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
# Down          1161             0                0          7051       3249     3     0
# NotSig        4215         20470            20470          8477      14966 20464 20470
# Up           15094             0                0          4942       2255     3     0

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
  geom_text_repel(data=subset(resFit_conditionvasc, abs(logFC)>1.5), aes(label=gene), col="black")+
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


pdf(file=paste0(dirQCplots,"volcanoPlot_vco.pdf"), width=8, height=7)
plot(volcanoPlot)
dev.off()

volcanoPlot <- ggplot(data=resFit_conditionvasc_mg, aes(x=logFC, y=-log10(adj.P.Val), col=signif))+
  geom_point(size=1, pch=15)+
  theme_bw()+
  scale_color_manual(name="DE genes",values=c("TRUE"="red","FALSE"="grey80"))+
  geom_hline(yintercept=-log10(0.1), linetype="dashed", col="grey30")+
  geom_vline(xintercept=c(-log2(1.5),-log2(2),-log2(3),log2(1.5), log2(2), log2(3), log2(4)), linetype="dashed", col="grey30")+
  geom_text_repel(data=subset(resFit_conditionvasc_mg, abs(logFC)>1.5), aes(label=gene), col="black")+
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


pdf(file=paste0(dirQCplots,"volcanoPlot_vico.pdf"), width=8, height=7)
plot(volcanoPlot)
dev.off()



MAplot <- ggplot(data=resFit_conditionvasc, aes(x=AveExpr, y=logFC, col=signif))+
  geom_point(size=1, pch=15, alpha=0.8)+
  theme_bw()+
  scale_color_manual(name="DE genes",values=c("TRUE"="red","FALSE"="grey80"))+
  geom_hline(yintercept=c(-log2(1.5),-log2(2),-log2(3),log2(1.5), log2(2), log2(3)), linetype="dashed", col="grey30")+
  geom_text_repel(data=subset(resFit_conditionvasc, abs(logFC)>1.5), aes(label=gene), col="black")+
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


pdf(file=paste0(dirQCplots,"MAPlot_proportionFit_vco.pdf"),  width=8, height=7)
plot(MAplot)
dev.off()

MAplot <- ggplot(data=resFit_conditionvasc_mg, aes(x=AveExpr, y=logFC, col=signif))+
  geom_point(size=1, pch=15, alpha=0.8)+
  theme_bw()+
  scale_color_manual(name="DE genes",values=c("TRUE"="red","FALSE"="grey80"))+
  geom_hline(yintercept=c(-log2(1.5),-log2(2),-log2(3),log2(1.5), log2(2), log2(3)), linetype="dashed", col="grey30")+
  geom_text_repel(data=subset(resFit_conditionvasc, abs(logFC)>1.5), aes(label=gene), col="black")+
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


pdf(file=paste0(dirQCplots,"MAPlot_proportionFit_vico.pdf"),  width=8, height=7)
plot(MAplot)
dev.off()

#######
#######

### GO enrichment - upregulated and downregulated (for the disease and the two treatment)

geneUniverse <- data.frame(symbol=rownames(se.filt))

library(org.Hs.eg.db)

tabCorr <- AnnotationDbi::select(org.Hs.eg.db, geneUniverse$symbol, "ENTREZID", "SYMBOL")
tabCorr <- tabCorr[!is.na(tabCorr$ENTREZID),]
tabCorr <- tabCorr[!duplicated(tabCorr$SYMBOL),]

geneUniverse$entrezid <- tabCorr[match(geneUniverse$symbol,tabCorr$SYMBOL),]$ENTREZID
geneUniverse <- subset(geneUniverse, !is.na(entrezid))


coefficients_to_runGO <- c("conditionvasc","conditionvasc_mg","donorAlz43cl2","donorPLCG1")


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



cov_goenrich <- sapply(coefficients_to_runGO, function(t) {
  
  print(t)
  
  resFit_coeff <- topTable(fit, coef=t, n=Inf)
  
  # upReg_list <- rownames(resFit_coeff[which(resFit_coeff$adj.P.Val<0.1 & resFit_coeff$logFC>0),])
  # downReg_list <- rownames(resFit_coeff[which(resFit_coeff$adj.P.Val<0.1 & resFit_coeff$logFC<0),])
  
  upReg_list <- rownames(resFit_coeff[which(resFit_coeff$logFC>1.5),])
  downReg_list <- rownames(resFit_coeff[which(resFit_coeff$logFC< -1.5),])
  
  list_up_down <- list(upReg_list,downReg_list)
  names(list_up_down) <- c("Upregulated","Downregulated")
  
  prop_go <- sapply(1:length(list_up_down), function(tt){
    
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
    
    prop_go$ranking <- 1:dim(prop_go)[1]
    #allres_intersect <- rbind(upreg_go[1:10,], downreg_go[1:10,])
    prop_go$short <- prop_go$Term
    mask_long<- nchar(prop_go$Term)>30
    
    return(prop_go)
  }, simplify=F)
  
  return(prop_go)
  
}, simplify=F)

cov_goenrich2 <- unlist(cov_goenrich, recursive = FALSE)

cov_goenrich_high <- cov_goenrich2
cov_goenrich_high <- sapply(cov_goenrich_high, function(x) x[1:10,], simplify=F)

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
  facet_wrap(
    ~label,
    nrow = 2,
    scales = "free_y",
    labeller = as_labeller(function(x) {
      x <- gsub("conditionvasc_mg", "<b>viCO</b>", x)
      x <- gsub("conditionvasc", "<b>vCO</b>", x)
      x <- gsub("donorAlz43cl2", "<b>CL2</b>", x)
      x <- gsub("donorPLCG1", "<b>CL3</b>", x)
      x
    })
    )+
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
  ggtitle("GO enrichment per covariate")

pdf(file=paste0(dirQCplots,"GOenrich_DE.pdf"), width=14.5, height=6)
plot(gostats_plot)
dev.off()





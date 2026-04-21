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
library(SeuratDisk)
library(tidyverse)
library(patchwork)


pathToOutput <- "/scratch/project_2012655/pau/saved/scanpy/output/"
branchTag <- "glutamatergicN"

alldonors <- dir(pathToOutput)[grepl("Tradeseq.glutamatergicN.", dir(pathToOutput))]
alldonors <- gsub(".tsv","",gsub("Tradeseq.glutamatergicN.","",alldonors))

rootDir <- "/scratch/project_2012655/pau/saved/scanpy/figures/"

sapply(alldonors, function(x){
  
  print(x)
  CombinedDF <- read.table(paste0(pathToOutput,"Tradeseq.glutamatergicN.",x,".tsv"), header=T)
  
  ## Top-upregulated and down-regulated genes
  topbot <- read.table(paste0(pathToOutput,"topbot_fitted.glutamatergicN.",x,".tsv"), header=T)
  colnames(topbot) <- gsub("fittedCounts\\.","", colnames(topbot))
  topbot <- as.data.frame(topbot)
  topbot$bucket <- paste0("b",1:dim(topbot)[1])
  topBotGenes <- as.data.frame(topbot %>% pivot_longer(-c("bucket"), names_to="genes", values_to="fittedExpression"))
  topBotGenes$genes <- gsub("\\.","-",topBotGenes$genes)
  topBotGenes$direction <- ifelse( CombinedDF[match(unique(topBotGenes$genes), CombinedDF$Gene),]$logFC>0, "Upregulated", "Downregulated")
  stopifnot(unique(topBotGenes$genes) %in% CombinedDF$Gene)
  topBotGenes$bucket <- factor(topBotGenes$bucket, levels= unique(topBotGenes$bucket))
  
  
  ## Cortical genes
  cortical <- read.table(paste0(pathToOutput,"cortical_relevant_fitted.glutamatergicN.",x,".tsv"), header=T)
  colnames(cortical) <- gsub("fittedCounts\\.","", colnames(cortical))
  cortical <- as.data.frame(cortical)
  cortical$bucket <- paste0("b",1:dim(cortical)[1])
  corticalGenes <- as.data.frame(cortical %>% pivot_longer(-c("bucket"), names_to="genes", values_to="fittedExpression"))
  corticalGenes$genes <- gsub("\\.","-",corticalGenes$genes)
  corticalGenes$direction <- ifelse( CombinedDF[match(unique(corticalGenes$genes), CombinedDF$Gene),]$logFC>0, "Upregulated", "Downregulated")
  corticalGenes$direction[is.na(corticalGenes$direction)] <- "No-DE"
  corticalGenes$bucket <- factor(corticalGenes$bucket, levels= unique(corticalGenes$bucket))
  
  ## PC1-relevant
  
  pc1 <- read.table(paste0(pathToOutput,"pc1_relevant_fitted.glutamatergicN.",x,".tsv"), header=T)
  colnames(pc1) <- gsub("fittedCounts\\.","", colnames(pc1))
  pc1 <- as.data.frame(pc1)
  pc1$bucket <- paste0("b",1:dim(pc1)[1])
  pc1Genes <- as.data.frame(pc1 %>% pivot_longer(-c("bucket"), names_to="genes", values_to="fittedExpression"))
  pc1Genes$genes <- gsub("\\.","-",pc1Genes$genes)
  pc1Genes$direction <- ifelse( CombinedDF[match(unique(pc1Genes$genes), CombinedDF$Gene),]$logFC>0, "Upregulated", "Downregulated")
  pc1Genes$direction[is.na(pc1Genes$direction)] <- "No-DE"
  pc1Genes$bucket <- factor(pc1Genes$bucket, levels= unique(pc1Genes$bucket))
  
  
  ## PC2-relevant
  
  pc2 <- read.table(paste0(pathToOutput,"pc2_relevant_fitted.glutamatergicN.",x,".tsv"), header=T)
  colnames(pc2) <- gsub("fittedCounts\\.","", colnames(pc2))
  pc2 <- as.data.frame(pc2)
  pc2$bucket <- paste0("b",1:dim(pc2)[1])
  pc2Genes <- as.data.frame(pc2 %>% pivot_longer(-c("bucket"), names_to="genes", values_to="fittedExpression"))
  pc2Genes$genes <- gsub("\\.","-",pc2Genes$genes)
  pc2Genes$direction <- ifelse( CombinedDF[match(unique(pc2Genes$genes), CombinedDF$Gene),]$logFC>0, "Upregulated", "Downregulated")
  pc2Genes$direction[is.na(pc2Genes$direction)] <- "No-DE"
  pc2Genes$bucket <- factor(pc2Genes$bucket, levels= unique(pc2Genes$bucket))
  
  
  
  all_genes <- list(
    unique(topBotGenes$genes),
    unique(corticalGenes$genes),
    unique(pc1Genes$genes),
    unique(pc2Genes$genes)
  )
  max_genes <- max(sapply(all_genes, length))  # maximum number of facets
  
  # For each dataset, add dummy genes if needed
  add_dummy_facets <- function(df, max_facets) {
    current_genes <- unique(df$genes)
    n_missing <- max_facets - length(current_genes)
    if(n_missing > 0){
      dummy <- data.frame(
        bucket = factor(rep("b1", n_missing)),  # dummy bucket
        genes = paste0("dummy_", seq_len(n_missing)),
        fittedExpression = NA,
        direction = "No-DE"
      )
      df <- rbind(df, dummy)
    }
    return(df)
  }
  
  topBotGenes <- add_dummy_facets(topBotGenes, max_genes)
  corticalGenes <- add_dummy_facets(corticalGenes, max_genes)
  pc1Genes <- add_dummy_facets(pc1Genes, max_genes)
  pc2Genes <- add_dummy_facets(pc2Genes, max_genes)
  
  # Custom labeller: blank for dummy facets, keep real gene names
  custom_labeller <- as_labeller(function(x){
    # x is a character vector of facet values
    sapply(x, function(val){
      if (grepl("^dummy_", val)) "" else val
    })
  })
  
  ### Plots
  
  topBotGenes$genes <- factor(topBotGenes$genes, levels=unique(topBotGenes$genes))
  corticalGenes$genes <- factor(corticalGenes$genes, levels=unique(corticalGenes$genes))
  pc1Genes$genes <- factor(pc1Genes$genes, levels=unique(pc1Genes$genes))
  pc2Genes$genes <- factor(pc2Genes$genes, levels=unique(pc2Genes$genes))
  
  
  topBotPlot <- ggplot(data = topBotGenes, aes(x = bucket, y = fittedExpression, fill=direction, col=direction)) +
    geom_bar(stat = "identity") +
    facet_wrap(~genes, scales = "free", nrow = 1, labeller = custom_labeller) +
    theme_void() +  # remove everything
    theme(
      strip.background = element_blank(),       # remove box around facet title
      strip.text = element_text(size = 8),     # keep facet title
      panel.spacing = unit(1, "lines"),         # optional spacing between panels
      plot.margin = margin(10, 10, 10, 10)      # optional outer margin
    )+
    scale_fill_manual(name="tradeSeq (SvsE)", values=c("Upregulated"="#E69F00","Downregulated"="#56B4E9", "No-DE"="grey75"))+
    scale_color_manual(name="tradeSeq (SvsE)", values=c("Upregulated"="#E69F00","Downregulated"="#56B4E9", "No-DE"="grey75"))+
    theme(legend.position="top",
          legend.title=element_text(size=8),
          legend.text=element_text(size=6),
          legend.key.size = unit(0.3, "cm"))+
    ggtitle("Top and Bot DE-Trajectory")
  
  
  corticalPlot <- ggplot(data = corticalGenes, aes(x = bucket, y = fittedExpression, fill=direction, col=direction)) +
    geom_bar(stat = "identity") +
    facet_wrap(~genes, scales = "free", nrow = 1, labeller = custom_labeller) +
    theme_void() +  # remove everything
    theme(
      strip.background = element_blank(),       # remove box around facet title
      strip.text = element_text(size = 8),     # keep facet title
      panel.spacing = unit(1, "lines"),         # optional spacing between panels
      plot.margin = margin(10, 10, 10, 10)      # optional outer margin
    )+
    scale_fill_manual(name="tradeSeq (SvsE)", values=c("Upregulated"="#E69F00","Downregulated"="#56B4E9", "No-DE"="grey75"))+
    scale_color_manual(name="tradeSeq (SvsE)", values=c("Upregulated"="#E69F00","Downregulated"="#56B4E9", "No-DE"="grey75"))+
    theme(legend.position="top",
          legend.title=element_text(size=8),
          legend.text=element_text(size=6),
          legend.key.size = unit(0.3, "cm"))+
    ggtitle("Cortical-relevant Genes")
  

  pc1Plot <- ggplot(data = pc1Genes, aes(x = bucket, y = fittedExpression, fill=direction, col=direction)) +
    geom_bar(stat = "identity") +
    facet_wrap(~genes, scales = "free", nrow = 1, labeller = custom_labeller) +
    theme_void() +  # remove everything
    theme(
      strip.background = element_blank(),       # remove box around facet title
      strip.text = element_text(size = 8),     # keep facet title
      panel.spacing = unit(1, "lines"),         # optional spacing between panels
      plot.margin = margin(10, 10, 10, 10)      # optional outer margin
    )+
    scale_fill_manual(name="tradeSeq (SvsE)", values=c("Upregulated"="#E69F00","Downregulated"="#56B4E9", "No-DE"="grey75"))+
    scale_color_manual(name="tradeSeq (SvsE)", values=c("Upregulated"="#E69F00","Downregulated"="#56B4E9", "No-DE"="grey75"))+
    theme(legend.position="top",
          legend.title=element_text(size=8),
          legend.text=element_text(size=6),
          legend.key.size = unit(0.3, "cm"))+
    ggtitle("PC1-Top/Bot Genes")
  
  

  pc2Plot <- ggplot(data = pc2Genes, aes(x = bucket, y = fittedExpression, fill=direction, col=direction)) +
    geom_bar(stat = "identity") +
    facet_wrap(~genes, scales = "free", nrow = 1, labeller = custom_labeller) +
    theme_void() +  # remove everything
    theme(
      strip.background = element_blank(),       # remove box around facet title
      strip.text = element_text(size = 8),     # keep facet title
      panel.spacing = unit(1, "lines"),         # optional spacing between panels
      plot.margin = margin(10, 10, 10, 10)      # optional outer margin
    )+
    scale_fill_manual(name="tradeSeq (SvsE)", values=c("Upregulated"="#E69F00","Downregulated"="#56B4E9", "No-DE"="grey75"))+
    scale_color_manual(name="tradeSeq (SvsE)", values=c("Upregulated"="#E69F00","Downregulated"="#56B4E9", "No-DE"="grey75"))+
    theme(legend.position="top",
          legend.title=element_text(size=8),
          legend.text=element_text(size=6),
          legend.key.size = unit(0.3, "cm"))+
    ggtitle("PC2-Top/Bot Genes")
  
  
 
  combinedPlot <- wrap_plots(topBotPlot, corticalPlot, pc1Plot, pc2Plot, nrow = 4) +
    plot_annotation(
      title = paste0("Gene expression patterns across pseudotime bins (",x,")"),
      theme = theme(
        plot.title = element_text(
          size = 14,
          face = "bold",
        )
      )
    )
  
  pdf(file=paste0(rootDir,"pseudotime_gene_patterns_",x,".pdf"), width=10, height = 6)
  plot(combinedPlot)
  dev.off()
  
  
})















library(tradeSeq)
library(ggplot2)
library(tidyverse)
library(scales)

set.seed(8)

pathToDir <- "saved/scanpy/"

w <- as.matrix(read.csv(paste0(pathToDir,"cellWeights_neurons2.csv"), header = FALSE))
dpt <- as.matrix(read.csv(paste0(pathToDir,"pseudotime_neurons2.csv"), header = FALSE))
cMatrix <- as.matrix(read.csv(paste0(pathToDir,"counts_neurons2.csv"), header = FALSE, row.names = 1))

gamObj <- fitGAM(cMatrix, verbose = TRUE, pseudotime = dpt, cellWeights = w, nknots = 6, sce=FALSE)
names(gamObj) <- rownames(cMatrix)

branchTag <- "glutamatergicN"

saveRDS(gamObj, file=paste0(pathToDir,"fitGAM_result_allgenes.",branchTag,".rds"))

dptseq <- seq(min(dpt), max(dpt), length.out = 5)

startRes <- startVsEndTest(gamObj, pseudotimeValues = c(min(dpt)+.01,max(dpt)-.01))

startResFilt <- startRes[startRes$pvalue <=0.05 & abs(startRes$logFClineage1)>=2,]
startResFilt$Gene <- rownames(startResFilt)
startResFilt$test <- "startVsEndTest"

nGenes <- 3

# Top-decreasing (positive to negative)
startResGenesPositive <- startResFilt[startResFilt$logFClineage1<0,]
startResGenesPositive <- startResGenesPositive[order(startResGenesPositive$pvalue,
                                                     startResGenesPositive$logFClineage1), ]
startResGenesPositive_plot <- startResGenesPositive[1:nGenes,]$Gene
startResGenesPositive$PatternType <- "decreasing"
colnames(startResGenesPositive)[colnames(startResGenesPositive)=="logFClineage1"] <- "logFC"

## Top-increasing (negative to positive)
startResGenesNegatives <- startResFilt[startResFilt$logFClineage1>0,]
startResGenesNegatives <- startResGenesNegatives[order(startResGenesNegatives$pvalue,
                                                       -startResGenesNegatives$logFClineage1), ]
startResGenesNegatives_plot <- startResGenesNegatives[1:nGenes,]$Gene
startResGenesNegatives$PatternType <- "increasing"
colnames(startResGenesNegatives)[colnames(startResGenesNegatives)=="logFClineage1"] <- "logFC"

rownames(startResGenesPositive) <- NULL
rownames(startResGenesNegatives) <- NULL

CombinedDF <- rbind(startResGenesPositive, startResGenesNegatives)
CombinedDF$Branch <- branchTag

write.table(CombinedDF, file=paste0(pathToDir,"Tradeseq.",branchTag,".tsv"), row.names = F, col.names = T, sep="\t", quote=F)

CombinedList <- c(startResGenesPositive_plot,startResGenesNegatives_plot)

AlternativeList <- c("SOX2","PAX6","HES1","EOMES","DCX","NEUROD1","STMN2","MAP2","SLC17A7","TBR1","CUX1")
pc1_relevant <- c("STMN2","MAPT","DCX","NRXN1","MAP2","SOX2","GLI3","VIM","HES1","CREB5")
pc2_relevant <- c("NFIB","NFIA","SOX5","PTPRZ1","SLC4A10","ERBB4","ZIC1","NR2F1","PAX3","FIGN")

.getPredictRangeDf <- function(dm, lineageId, conditionId = NULL, nPoints = 50){
  vars <- dm[1, ]
  if ("y" %in% colnames(vars)) {
    vars <- vars[!colnames(vars) %in% "y"]
    off <- 1
  } else {
    off <- 0
  }
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # duplicate to nPoints
  vars <- rbind(vars, vars[rep(1, nPoints - 1), ])
  rownames(vars) <- 1:length(rownames(vars))
  
  # set range of pseudotime for lineage of interest
  if (is.null(conditionId)) {
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId))  
  } else {
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, conditionId))
  }
  if (length(lineageIds) == 1){
    lineageData <- dm[dm[, lineageIds + off] == 1,
                      paste0("t", lineageId)]
  } else {
    lineageData <- dm[rowSums(dm[, lineageIds + off]) == 1,
                      paste0("t", lineageId)]
  }
  
  # make sure lineage starts at zero
  if(min(lineageData) / max(lineageData) < .01) {
    lineageData[which.min(lineageData)] <- 0
  }
  
  vars[, lineageIds] <- 1 / length(lineageIds)
  # set lineage
  vars[, paste0("t", lineageId)] <- seq(min(lineageData),
                                        max(lineageData),
                                        length = nPoints)
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                       pattern = "offset")])
  return(vars)
}



datalist <- list()
for (g in CombinedList){

  localModel <- gamObj[[g]]
  data <- localModel$model
  y <- data$y
  nCurves <- length(localModel$smooth)

  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(localModel$model, jj, nPoints = 100)
    yhat <- predict(localModel, newdata = df, type = "response")

    Newframe <- data.frame("fittedCounts" = yhat)
    colnames(Newframe) <- c(paste0("fittedCounts.",g))

    datalist[[paste(g,jj)]] <- Newframe

  }

}


datalist2 <- list()
for (g in AlternativeList){
  
  localModel <- gamObj[[g]]    
  data <- localModel$model    
  y <- data$y
  nCurves <- length(localModel$smooth)
  
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(localModel$model, jj, nPoints = 100)
    yhat <- predict(localModel, newdata = df, type = "response")
    
    Newframe <- data.frame("fittedCounts" = yhat)
    colnames(Newframe) <- c(paste0("fittedCounts.",g))
    
    datalist2[[paste(g,jj)]] <- Newframe
    
  }
  
}


datalist3 <- list()
for (g in pc1_relevant){
  
  localModel <- gamObj[[g]]    
  data <- localModel$model    
  y <- data$y
  nCurves <- length(localModel$smooth)
  
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(localModel$model, jj, nPoints = 100)
    yhat <- predict(localModel, newdata = df, type = "response")
    
    Newframe <- data.frame("fittedCounts" = yhat)
    colnames(Newframe) <- c(paste0("fittedCounts.",g))
    
    datalist3[[paste(g,jj)]] <- Newframe
    
  }
  
}


datalist4 <- list()
for (g in pc2_relevant){
  
  localModel <- gamObj[[g]]    
  data <- localModel$model    
  y <- data$y
  nCurves <- length(localModel$smooth)
  
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(localModel$model, jj, nPoints = 100)
    yhat <- predict(localModel, newdata = df, type = "response")
    
    Newframe <- data.frame("fittedCounts" = yhat)
    colnames(Newframe) <- c(paste0("fittedCounts.",g))
    
    datalist4[[paste(g,jj)]] <- Newframe
    
  }
  
}


topbot_fitted =  do.call(cbind, datalist)
cortical_relevant_fitted <- do.call(cbind, datalist2)
pc1_relevant_fitted <- do.call(cbind, datalist3)
pc2_relevant_fitted <- do.call(cbind, datalist4)


write.table(topbot_fitted, file=paste0(pathToDir,"topbot_fitted.",branchTag,".tsv"), row.names = T, col.names = T, sep="\t", quote=F)
write.table(cortical_relevant_fitted, file=paste0(pathToDir,"cortical_relevant_fitted.",branchTag,".tsv"), row.names = T, col.names = T, sep="\t", quote=F)
write.table(pc1_relevant_fitted, file=paste0(pathToDir,"pc1_relevant_fitted.",branchTag,".tsv"), row.names = T, col.names = T, sep="\t", quote=F)
write.table(pc2_relevant_fitted, file=paste0(pathToDir,"pc2_relevant_fitted.",branchTag,".tsv"), row.names = T, col.names = T, sep="\t", quote=F)


#########
#########

CombinedDF <- read.table(paste0(pathToDir,"Tradeseq.",branchTag,".tsv"), header = T)

## Top-upregulated and down-regulated genes
topbot <- read.table(paste0(pathToDir,"topbot_fitted.",branchTag,".tsv"))
colnames(topbot) <- gsub("fittedCounts\\.","", colnames(topbot))
topbot <- as.data.frame(topbot)
topbot$bucket <- paste0("b",1:dim(topbot)[1])
topBotGenes <- as.data.frame(topbot %>% pivot_longer(-c("bucket"), names_to="genes", values_to="fittedExpression"))
topBotGenes$genes <- gsub("\\.","-",topBotGenes$genes)

upGenes <- c("S1PR1","H3C2","H2BC9")
downGenes <- c("LINC01551","ENSG00000257522","NEUROD6")

topBotGenes$direction <- ifelse(CombinedDF[match(topBotGenes$genes, CombinedDF$Gene),]$logFC>0, "Upregulated", "Downregulated")
stopifnot(unique(topBotGenes$genes) %in% CombinedDF$Gene)
topBotGenes$bucket <- factor(topBotGenes$bucket, levels= unique(topBotGenes$bucket))


## Cortical genes
cortical <- read.table(paste0(pathToDir,"cortical_relevant_fitted.",branchTag,".tsv"), header=T)
colnames(cortical) <- gsub("fittedCounts\\.","", colnames(cortical))
cortical <- as.data.frame(cortical)
cortical$bucket <- paste0("b",1:dim(cortical)[1])
corticalGenes <- as.data.frame(cortical %>% pivot_longer(-c("bucket"), names_to="genes", values_to="fittedExpression"))
corticalGenes$genes <- gsub("\\.","-",corticalGenes$genes)
corticalGenes$direction <- ifelse( CombinedDF[match(unique(corticalGenes$genes), CombinedDF$Gene),]$logFC>0, "Upregulated", "Downregulated")
corticalGenes$direction[is.na(corticalGenes$direction)] <- "No-DE"
corticalGenes$bucket <- factor(corticalGenes$bucket, levels= unique(corticalGenes$bucket))

## PC1-relevant

pc1 <- read.table(paste0(pathToDir,"pc1_relevant_fitted.",branchTag,".tsv"), header=T)
colnames(pc1) <- gsub("fittedCounts\\.","", colnames(pc1))
pc1 <- as.data.frame(pc1)
pc1$bucket <- paste0("b",1:dim(pc1)[1])
pc1Genes <- as.data.frame(pc1 %>% pivot_longer(-c("bucket"), names_to="genes", values_to="fittedExpression"))
pc1Genes$genes <- gsub("\\.","-",pc1Genes$genes)
pc1Genes$direction <- ifelse( CombinedDF[match(unique(pc1Genes$genes), CombinedDF$Gene),]$logFC>0, "Upregulated", "Downregulated")
pc1Genes$direction[is.na(pc1Genes$direction)] <- "No-DE"
pc1Genes$bucket <- factor(pc1Genes$bucket, levels= unique(pc1Genes$bucket))


## PC2-relevant

pc2 <- read.table(paste0(pathToDir,"pc2_relevant_fitted.",branchTag,".tsv"), header=T)
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
    title = paste0("Gene expression patterns across pseudotime bins (All cells, Lines&conditions mixed)"),
    theme = theme(
      plot.title = element_text(
        size = 14,
        face = "bold",
      )
    )
  )

pdf(file=paste0(rootDir,"pseudotime_gene_patterns_allCells.pdf"), width=10, height = 6)
plot(combinedPlot)
dev.off()



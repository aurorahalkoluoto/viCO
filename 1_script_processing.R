


library(Seurat)
library(ggplot2)
library(scales)

list_10xfiles <- read.table("/scratch/project_2012655/all10xlibs.txt")$V1

processing <- function(list_10xfiles, return.obj=T, mito=F){
  
  allobj <- sapply(list_10xfiles, function(x){
    print(paste0("Processing sample ",x))
    tenXrun <- Read10X(data.dir = x)[[1]]
    
    sample <- sapply(strsplit(x,"/"), function(x) x[5])
    
    tenXrun_seurat <- CreateSeuratObject(counts = tenXrun, project =sample)
    print(dim(tenXrun_seurat))
    tenXrun_seurat@meta.data$donorId <-  sample
    rm(tenXrun)
    gc()
    
    ## NPC were frozen at day 14 at Sanger.
    ## Then, they were
    
    limsNFeatures <- c(2500,7500)
    limsMito <- 10
    
    idDef <- unique(tenXrun_seurat@meta.data$donorId)
    tenXrun_seurat[["percent.mt"]] <- PercentageFeatureSet(tenXrun_seurat, pattern = "^MT-")
    
    plot1 <- FeatureScatter(tenXrun_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")+
      ggtitle(idDef)+theme(plot.title=element_text(hjust=0.5, size=11),
                           legend.position="none")
    plot2 <- FeatureScatter(tenXrun_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
      ggtitle(idDef)+
      theme(plot.title=element_text(hjust=0.5, size=11),
            legend.position="none")
    qcplot1 <- plot1 + plot2
    
    tenXrun_seurat[["percent.ribo"]]<- PercentageFeatureSet(tenXrun_seurat, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
    ##QC-ribosomal content
    stopifnot(all(mean(tenXrun_seurat@meta.data$percent.ribo)>10 | mean(tenXrun_seurat@meta.data$percent.ribo)<20))
    
    dirQCplots <- paste0("/scratch/project_2012655/pau/QCplots")
    pdf(file=paste0(dirQCplots,"QCplot_",idDef, ".pdf"))
    plot(qcplot1)
    dev.off()
    
    tmpTabQCFeature <- data.frame(nFeature_RNA=unname(sort(tenXrun_seurat$nFeature_RNA)),
                                  index=1:length(sort(tenXrun_seurat$nFeature_RNA)),
                                  sample=sapply(strsplit(x, "/"), function(y) y[5]))
    
    tmpTabQCMito <- data.frame(percent.mt=unname(sort(tenXrun_seurat$percent.mt)),
                               index=1:length(sort(tenXrun_seurat$percent.mt)),
                               sample=sapply(strsplit(x, "/"), function(y) y[5]))
    
    
    tenXrun_seurat <- subset(tenXrun_seurat, subset = nFeature_RNA > limsNFeatures[1] & nFeature_RNA < limsNFeatures[2] & percent.mt < limsMito)
    gc()
    
    if (return.obj==TRUE){
      print(dim(tenXrun_seurat))
      return(tenXrun_seurat)
      
    } else {
      
      if (mito==TRUE){
        
        return(tmpTabQCMito)
        
      } else {
        
        return(tmpTabQCFeature)
        
      }
    }
    
  }, simplify=F)
  
  
  return(allobj)
  
  
}

allobj <- processing(list_10xfiles, return.obj=T, mito=F)
names(allobj) <- sapply(strsplit(names(allobj),"/"), function(x) x[length(x)])



saveRDS(allobj, file="/scratch/project_2012655/pau/allLibraries_notMerged_list.RDS")


nCountQC <- processing(list_10xfiles, return.obj=F, mito=F)
mitoQC <- processing(list_10xfiles, return.obj=F, mito=T)

nCountQC <- do.call("rbind", nCountQC)
rownames(nCountQC) <- NULL
mitoQC <- do.call("rbind", mitoQC)
rownames(mitoQC) <- NULL


saveRDS(nCountQC, file="/scratch/project_2012655/pau/qc_nFeat.RDS")
saveRDS(mitoQC, file="/scratch/project_2012655/pau/qc_mito.RDS")


featuresPlot <- ggplot(nCountQC, aes(x=index, y=nFeature_RNA))+
  geom_line()+
  theme_bw()+
  theme(legend.position="top",
        plot.title=element_text(hjust=0.5, face="bold", size=13))+
  theme(legend.title = element_text( size=8), legend.text=element_text(size=8))+
  scale_color_viridis_d(name="")+
  ggtitle("Mean number of expressed genes per cell")+
  geom_hline(yintercept=c(2500,7500), col="red", linetype="dashed")+
  scale_y_continuous(labels = comma, limits=c(0,12500),breaks=seq(0,12500,2500))+
  xlab("Cells ordered by number of genes expressed (in each 10x Sample)")+
  ylab("")+
  facet_wrap(~sample)


#mitoQC$sample2 <- gsub("OT_NPC23-1_","", mitoQC$sample)

mitoPlot <- ggplot(mitoQC, aes(x=index, y=percent.mt))+
  geom_line()+
  theme_bw()+
  theme(legend.position="top",
        plot.title=element_text(hjust=0.5, face="bold", size=13))+
  guides(col=guide_legend(nrow=4,byrow=TRUE))+
  theme(legend.title = element_text( size=8), legend.text=element_text(size=8))+
  scale_color_viridis_d(name="")+
  ggtitle("% of Mitochondrial Reads")+
  geom_hline(yintercept=c(10), col="red", linetype="dashed")+
  scale_y_continuous(labels = comma, limits=c(0,100),breaks=seq(0,100,10))+
  xlab("Cells ordered by % of mito reads (in each 10x Sample)")+
  ylab("")+
  facet_wrap(~sample)

dirQCplots <- "/scratch/project_2012655/pau/QCplots/"
pdf(file=paste0(dirQCplots,"QC_numGenesPerCell.pdf"), width=10, height = 8)
plot(featuresPlot)
dev.off()

pdf(file=paste0(dirQCplots,"QC_mitoPerc.pdf"), width=10, height = 8)
plot(mitoPlot)
dev.off()


## merge all 10x libraries (with raw counts, filtered by nFeatures and mito counts. Starting object for any downstream work)


mergedObj <- merge(allobj[[1]], y=do.call("c",allobj[2:length(allobj)]),
                   add.cell.ids = names(allobj),
                   project="Organoids")
gc()

mergedObj <- JoinLayers(mergedObj)

## genes expressed in less than 0.1% total cells are removed
selected_f <- names(which(rowSums(mergedObj@assays$RNA["counts"]>0)>0.005*dim(mergedObj)[2]))
mergedObj <- subset(mergedObj, features = selected_f)
saveRDS(mergedObj, file="/scratch/project_2012655/pau/Organoids_MergedRawCounts.RDS")





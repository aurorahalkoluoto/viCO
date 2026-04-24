# Sparse pseudobulk korrelaatio atlas vs organoid 
# Fig2D, Sup1B


library(Matrix)
library(Seurat)
library(dplyr)
library(batchelor)
library(pheatmap)
library(grid)


cortex_ref <- NormalizeData(ref)


common_genes <- intersect(rownames(cortex_ref), rownames(all))
atlas <- subset(cortex_ref, features = common_genes)
org <- subset(all, features = common_genes)
all_cells <- colnames(all)
org_cells <- colnames(org)
org$condition <- all$condition[match(org_cells, all_cells)]
gc()


seurat.list <- list(atlas = atlas, organoid = org)
seurat.list <- lapply(seurat.list, FindVariableFeatures, nfeatures = 2000)
hvg <- Reduce(intersect, lapply(seurat.list, VariableFeatures))
gc()

mnn.out <- fastMNN(
  GetAssayData(seurat.list$atlas, layer = "counts")[hvg, ],
  GetAssayData(seurat.list$org, layer = "counts")[hvg, ]
)
gc()

reconstructed_mat <- assay(mnn.out, "reconstructed")
reconstructed_mat <- as(reconstructed_mat, "dgCMatrix")
integrated <- CreateSeuratObject(counts = reconstructed_mat)


integrated$dataset <- c(rep("atlas", ncol(atlas)), rep("organoid", ncol(org)))
integrated$CellClass <- c(atlas$CellClass, org$sctype_pred)


DefaultAssay(integrated) <- "RNA"

atlas_pb <- AverageExpression(
  subset(integrated, dataset == "atlas"),
  group.by = "CellClass",
  assays = "RNA",
  slot = "counts"
)$RNA

org_pb <- AverageExpression(
  subset(integrated, dataset == "organoid"),
  group.by = "CellClass",
  assays = "RNA",
  slot = "counts"
)$RNA


common_pb_genes <- intersect(rownames(atlas_pb), rownames(org_pb))
atlas_mat <- atlas_pb[common_pb_genes, , drop = FALSE]
org_mat <- org_pb[common_pb_genes, , drop = FALSE]


corr_mat <- cor(as.matrix(atlas_mat), as.matrix(org_mat), method = "pearson")
corr_mat_t <- t(corr_mat)
rownames(corr_mat_t) <- colnames(org_mat)
colnames(corr_mat_t) <- colnames(atlas_mat)


pheatmap(
  corr_mat_t,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = " ",
  angle_col = 45,
  fontsize_row = 14,  
  fontsize_col = 14    
)
grid::grid.text(
  "Fetal Cortex Atlas (pcw 6.9), Braun et al. (2023)",
  y = unit(0.013, "npc"),
  gp = grid::gpar(fontsize = 18)
)


grid::grid.text(
  "Integrated Organoid Dataset",
  x = unit(0.97, "npc"),
  y = unit(0.425, "npc"),
  rot = 270,
  gp = grid::gpar(fontsize = 18)
)




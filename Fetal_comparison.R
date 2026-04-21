Sparse pseudobulk korrelaatio atlas vs organoid
# =====================================================

library(Matrix)
library(Seurat)
library(dplyr)
library(batchelor)
library(pheatmap)
library(grid)


cortex_ref <- NormalizeData(ref_6)
#varmista että all annotoitu

#all <- subset(all, subset = sampletype == "vasc_mg")

# ---------- 4. Geenien yhdenmukaistus HVG-pohjalta ----------
common_genes <- intersect(rownames(cortex_ref), rownames(all))
atlas <- subset(cortex_ref, features = common_genes)
org <- subset(all, features = common_genes)
all_cells <- colnames(all)
org_cells <- colnames(org)
org$condition <- all$condition[match(org_cells, all_cells)]
gc()

# Etsi HVG:t
seurat.list <- list(atlas = atlas, organoid = org)
seurat.list <- lapply(seurat.list, FindVariableFeatures, nfeatures = 2000)
hvg <- Reduce(intersect, lapply(seurat.list, VariableFeatures))
gc()
# ---------- 5. fastMNN-integraatio ----------
mnn.out <- fastMNN(
  GetAssayData(seurat.list$atlas, layer = "counts")[hvg, ],
  GetAssayData(seurat.list$org, layer = "counts")[hvg, ]
)
gc()

# Luo Seurat-objekti suoraan reconstructed-matriisista
reconstructed_mat <- assay(mnn.out, "reconstructed")
reconstructed_mat <- as(reconstructed_mat, "dgCMatrix")
integrated <- CreateSeuratObject(counts = reconstructed_mat)

# ---------- 6. Lisää metadata ----------
integrated$dataset <- c(rep("atlas", ncol(atlas)), rep("organoid", ncol(org)))
integrated$CellClass <- c(atlas$CellClass, org$sctype_pred)

# ---------- 7. Pseudobulkit ----------
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

# ---------- 8. Yhdenmukaista geenit pseudobulkissa ----------
common_pb_genes <- intersect(rownames(atlas_pb), rownames(org_pb))
atlas_mat <- atlas_pb[common_pb_genes, , drop = FALSE]
org_mat <- org_pb[common_pb_genes, , drop = FALSE]

# ---------- 9. Pearson-korrelaatio ----------
corr_mat <- cor(as.matrix(atlas_mat), as.matrix(org_mat), method = "pearson")
corr_mat_t <- t(corr_mat)
rownames(corr_mat_t) <- colnames(org_mat)
colnames(corr_mat_t) <- colnames(atlas_mat)

# ---------- 10. Heatmap ----------
pheatmap(
  corr_mat_t,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = " ",
  angle_col = 45,
  fontsize_row = 14,   # Y-akselin (rivit) tekstit isommaksi
  fontsize_col = 14    # X-akselin (sarakkeet) tekstit isommaksi
)
grid::grid.text(
  "Fetal Cortex Atlas (pcw 6.9), Braun et al. (2023)",
  y = unit(0.013, "npc"),
  gp = grid::gpar(fontsize = 18)
)

# Y-akselin otsikko (vasemmalle, pystysuunnassa)
grid::grid.text(
  "Integrated Organoid Dataset",
  x = unit(0.97, "npc"),
  y = unit(0.425, "npc"),
  rot = 270,
  gp = grid::gpar(fontsize = 18)
)

# ---------- 11. Tallenna korrelaatiot ----------
write.csv(corr_mat, "atlas_pcw6_vs_organoid_fastMNN_correlation_v14_05.csv")
cat("✅ Heatmap ja korrelaatiot valmiit.\n")



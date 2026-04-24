#Fig 2B

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(ggrepel)




markers_list <- list( 
  ChP_epithelial = c("TTR", "KRT18", "KCNJ13"), 
  Endothelial = c("PECAM1", "CD34", "KDR"), 
  Fibroblast = c("COL1A1", "DCN", "LUM"), 
  GABAergic_neuron = c("GAD1", "GAD2", "SLC32A1"), 
  Glutamatergic_neuron = c("SLC17A7", "GRIN2B", "GRIN1"), 
  Microglia = c("P2RY12", "CX3CR1", "AIF1"),
  Neuroblast_Immature = c("DCX", "NEUROD1", "TUBB3"),
  Radial_glial = c("PAX6", "NES", "FABP7", "SOX2", "SLC1A3", "EOMES"),
  Roof_plate_progenitors = c("LMX1A", "MSX1", "BMP7"),
  Schwann_precursor = c("SOX10", "ERBB3", "NGFR"),
  Pericyte = c("RGS5", "PDGFRB", "MCAM"), 
  Smooth_muscle = c("ACTA2", "MYH11", "TAGLN")
  
)


marker_vector <- unlist(markers_list, use.names = FALSE)
marker_levels <- unique(marker_vector)

all$sctype_pred <- factor(
  all$sctype_pred,
  levels = sort(unique(all$sctype_pred))
)
Idents(all) <- "sctype_pred"

dp <- DotPlot(all, features = marker_levels)
dp_data <- dp$data


dp_data$gene <- factor(
  dp_data$features.plot,
  levels = rev(marker_levels)   
)


n_genes <- nrow(dp_data)  

hlines <- seq(3.5, n_genes - 0.5, by = 3)


p <- ggplot(dp_data, aes(x = id, y = gene)) +
  geom_point(aes(size = pct.exp, colour = avg.exp)) +
  scale_colour_gradient(
    low = "blue", high = "red",
    limits = c(0, 5),
    oob = scales::squish
  ) +
  geom_hline(
    yintercept = hlines,
    colour = "black",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(colour = "black"),
    axis.title.y = element_text(size = 10)
  ) +
  labs(
    x = "Cell type",
    y = "Marker genes",
    size = "% expressed",
    colour = "Avg expression"
  )

p



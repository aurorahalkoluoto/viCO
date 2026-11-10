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


#P2RY12 analysis

library(Seurat)
library(ggplot2)
library(patchwork)

# Set organoid type order
all$condition <- recode(
  all$condition,
  "ctr" = "cCO",
  "vasc" = "vCO",
  "vasc_mg" = "viCO"
)

all$condition <- factor(
  all$condition,
  levels = c("cCO", "vCO", "viCO")
)

# Plot P2RY12 expression separately for each organoid type
# UMAP plots
p_D <- FeaturePlot(
  all,
  features = "P2RY12",
  split.by = "condition",
  cols = c("grey85", "blue"),
  order = TRUE,
  pt.size = 0.3,
  min.cutoff = "q05",
  max.cutoff = "q95",
  combine = TRUE
) &
  guides(color = "none", fill = "none") &
  theme_classic()
p_D

p_E <- DotPlot(
  all,
  features = "P2RY12",
  group.by = "condition",
  cols = c("grey85", "blue"),
  dot.scale = 8
) +
  scale_y_discrete(limits = rev(c("viCO", "vCO", "cCO"))) +
  labs(
    title = "P2RY12 expression by organoid type",
    x = "Features",
    y = "Identity",
    color = "Average Expression",
    size = "Percent Expressed"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_E



library(dplyr)
library(tidyr)
library(ggplot2)

# Fetch data
p2ry12_data <- FetchData(
  all,
  vars = c("P2RY12", "condition", "sctype_pred"),
  layer = "data"
)

# Calculate percentages
p2ry12_summary <- p2ry12_data %>%
  group_by(condition) %>%
  summarise(
    in_all_cells = mean(P2RY12 > 0) * 100,
    in_microglia = mean(P2RY12[sctype_pred == "Microglia"] > 0) * 100,
    .groups = "drop"
  )

# Long format for ggplot
plot_data <- p2ry12_summary %>%
  pivot_longer(
    cols = c(in_all_cells, in_microglia),
    names_to = "group",
    values_to = "percent"
  )

# Plot
p_F <- ggplot(
  plot_data,
  aes(x = condition, y = percent, fill = group)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = round(percent, 1)),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "in_all_cells" = "grey75",
      "in_microglia" = "blue"
    ),
    labels = c(
      "in_all_cells",
      "in_microglia"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "% P2RY12+ cells",
    x = NULL,
    y = "% positive cells",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 11)
  )

p_F





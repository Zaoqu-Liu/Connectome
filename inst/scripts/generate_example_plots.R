# Script to generate example plots for vignettes
# Run this script to regenerate all example figures

library(Connectome)
library(ggplot2)
library(igraph)

# Create mock connectome data for visualization
set.seed(42)

cell_types <- c("Epithelial", "Fibroblast", "Macrophage", "T_cell", "Endothelial")
ligands <- c("VEGFA", "IL6", "TNF", "CXCL12", "FGF2", "TGFB1", "WNT5A", "BMP4")
receptors <- c("KDR", "IL6R", "TNFRSF1A", "CXCR4", "FGFR1", "TGFBR1", "FZD5", "BMPR1A")
modes <- c("growth_factor", "cytokine", "cytokine", "chemokine", "growth_factor", "growth_factor", "WNT", "BMP")

# Generate mock connectome
n_edges <- 200
mock_connectome <- data.frame(
  source = sample(cell_types, n_edges, replace = TRUE),
  target = sample(cell_types, n_edges, replace = TRUE),
  ligand = sample(ligands, n_edges, replace = TRUE),
  receptor = sample(receptors, n_edges, replace = TRUE),
  mode = sample(modes, n_edges, replace = TRUE),
  ligand.expression = runif(n_edges, 0.5, 3),
  recept.expression = runif(n_edges, 0.5, 3),
  ligand.scale = rnorm(n_edges, 0, 1),
  recept.scale = rnorm(n_edges, 0, 1),
  percent.source = runif(n_edges, 0.1, 0.8),
  percent.target = runif(n_edges, 0.1, 0.8),
  weight_norm = runif(n_edges, 0.5, 5),
  weight_sc = runif(n_edges, -2, 2)
)
mock_connectome$pair <- paste(mock_connectome$ligand, mock_connectome$receptor, sep = "_")
mock_connectome$vector <- paste(mock_connectome$source, mock_connectome$target, sep = " - ")

# Filter for visualization
mock_filtered <- mock_connectome[mock_connectome$weight_norm > 2, ]

# Output directory
outdir <- "vignettes/figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 1. Network Plot
message("Generating Network Plot...")
png(file.path(outdir, "network_plot.png"), width = 800, height = 600, res = 100)
NetworkPlot(mock_filtered, weight.attribute = "weight_norm", title = "Cell-Cell Communication Network")
dev.off()

# 2. Circos Plot  
message("Generating Circos Plot...")
png(file.path(outdir, "circos_plot.png"), width = 800, height = 800, res = 100)
CircosPlot(mock_filtered, weight.attribute = "weight_norm")
dev.off()

# 3. Edge Dot Plot
message("Generating Edge Dot Plot...")
p <- EdgeDotPlot(mock_filtered)
ggsave(file.path(outdir, "edge_dotplot.png"), p, width = 10, height = 8, dpi = 100)

# 4. Centrality Plot
message("Generating Centrality Plot...")
p <- Centrality(mock_filtered, group.by = "mode")
ggsave(file.path(outdir, "centrality_plot.png"), p, width = 10, height = 6, dpi = 100)

# 5. Modal Dot Plot
message("Generating Modal Dot Plot...")
p <- ModalDotPlot(mock_filtered, weight.attribute = "weight_norm")
ggsave(file.path(outdir, "modal_dotplot.png"), p, width = 12, height = 6, dpi = 100)

# 6. Signal Scatter
message("Generating Signal Scatter...")
p <- SignalScatter(mock_filtered, features = c("VEGFA", "IL6"), weight.attribute = "weight_norm")
ggsave(file.path(outdir, "signal_scatter.png"), p, width = 8, height = 6, dpi = 100)

message("All plots generated successfully!")
message(paste("Output directory:", normalizePath(outdir)))

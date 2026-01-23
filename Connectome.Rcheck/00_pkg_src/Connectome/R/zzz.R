#' @importFrom stats ave reorder
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom grDevices colorRampPalette recordPlot
#' @importFrom graphics par
#' @importFrom ggplot2 ggplot aes geom_point geom_tile theme theme_bw theme_light
#'   element_blank element_text element_rect ggtitle xlab ylab guides guide_legend
#'   scale_colour_manual scale_color_gradientn scale_fill_gradient2
#'   scale_x_discrete scale_y_discrete coord_flip ylim
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid get_legend
#' @importFrom igraph graph_from_data_frame hub_score authority_score
#'   layout_in_circle layout_with_fr V V<- E E<- get.edge.attribute
#' @importFrom circlize circos.clear chordDiagram circos.trackPlotRegion
#'   get.cell.meta.data circos.text
#' @importFrom ComplexHeatmap Legend draw
#' @importFrom grid gpar unit
#' @importFrom scales hue_pal
#' @importFrom dplyr %>% select
#' @importFrom gtools foldchange foldchange2logratio
#' @importFrom Seurat Idents Idents<- WhichCells SplitObject FindMarkers
#'   AverageExpression FindAllMarkers GetAssayData FetchData CreateSeuratObject
#'   CreateAssayObject DefaultAssay DefaultAssay<- NormalizeData ScaleData Assays
#' @importFrom Matrix rowSums
#' @importFrom data.table data.table as.data.table setorder rbindlist setcolorder
#' @keywords internal
"_PACKAGE"

# Global variable declarations to avoid R CMD check NOTEs
# These are column names used in data.frame/data.table operations
utils::globalVariables(c(

  # Column names from connectome data structure
  "source", "target", "ligand", "receptor", "pair", "mode", "vector", "edge",
  "ligand.expression", "recept.expression", "ligand.scale", "recept.scale",
  "percent.source", "percent.target", "weight_norm", "weight_sc",
  "DOR.source", "DOR.target", "p_val_adj", "p_val_adj.lig", "p_val_adj.rec",
  

  # Column names from differential connectome
  "ligand.norm.lfc", "recept.norm.lfc", "score",
  "pct.source.1", "pct.source.2", "pct.target.1", "pct.target.2",
  
  # data.table internal variables
  "id", ".",
  
  # Variables used in network analysis
  "group", "name", "cells", "hub.score", "auth.score",
  "wt.source", "wt.sink",
  
  # Future package variables
  "plan", "multisession",
  
  # ggplot aesthetic variables  
  "..count.."
))

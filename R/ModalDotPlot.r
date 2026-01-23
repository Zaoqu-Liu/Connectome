#' ModalDotPlot
#'
#' This function takes a connectomic edgelist and creates a source and a sink mode- and cell- organized dot plot. The y-axis is the discrete variable 'mode',
#' and the x-axis is the sum of the weights of all edges for each mode made by each cell. Points are organized by cell type, with the size of the point
#' correlating to the Kleinberg hub score (for source graph) and Kleinberg authority score (for sink graph). Network filtration is performed prior to network centrality calculations.
#'
#' @param connectome A connectomic edgelist
#' @param cols.use Desired colors for cell types, alphabetized. Defaults to standard ggplot colors.
#' @param weight.attribute Column to use to define edgeweights for network analysis.  'weight_sc' or 'weight_norm'. Defaults to 'weight_sc'
#' @param min.z Minimum z-score for ligand and receptor.
#' @param normalize Default TRUE. Scales each mode to have equivalent x-axes.
#' @param ... Arguments passed to FilterConnectome
#' @return A ggplot object showing modal centrality analysis
#' @export

ModalDotPlot <- function(connectome,
                         cols.use = NULL,
                         weight.attribute = "weight_sc",
                         min.z = NULL,
                         normalize = TRUE, ...) {
  # Check required packages
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required but not installed.")
  }

  if (weight.attribute == "weight_sc" && is.null(min.z)) {
    message("\nWeight attribute is 'weight_sc', recommend also setting min.z = 0 to avoid negative ligand and receptor scores")
  }

  # Store
  master <- connectome

  # Filter as demanded (remove NAs at minimum)
  master_sub <- FilterConnectome(master, min.z = min.z, remove.na = TRUE, ...)

  # Set up to plot ModalDotPlot
  modes <- as.character(unique(master_sub$mode))
  # Use unique values from data, not factor levels (more robust)
  cells <- as.character(unique(c(master_sub$source, master_sub$target)))

  # Process data
  n_modes <- length(modes)
  n_cells <- length(cells)
  result_list <- vector("list", n_modes * n_cells)
  idx <- 1

  for (i in seq_along(modes)) {
    temp <- master_sub[master_sub$mode == modes[[i]], ]
    
    if (nrow(temp) == 0) next
    
    net <- igraph::graph_from_data_frame(temp, directed = TRUE)
    hub <- igraph::hub_score(net, weights = temp[, weight.attribute], scale = TRUE)$vector
    auth <- igraph::authority_score(net, weights = temp[, weight.attribute], scale = TRUE)$vector
    total.edgeweight <- sum(temp[, weight.attribute])

    for (j in seq_along(cells)) {
      wt.source <- sum(temp[temp$source == cells[[j]], weight.attribute])
      wt.sink <- sum(temp[temp$target == cells[[j]], weight.attribute])
      
      if (normalize && total.edgeweight > 0) {
        wt.source <- wt.source / total.edgeweight
        wt.sink <- wt.sink / total.edgeweight
      }

      result_list[[idx]] <- data.frame(
        mode = modes[[i]], 
        cells = cells[[j]], 
        hub.score = hub[cells[[j]]], 
        auth.score = auth[cells[[j]]], 
        wt.source = wt.source, 
        wt.sink = wt.sink,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  
  # Combine all results at once
  df <- do.call(rbind, result_list[seq_len(idx - 1)])
  
  # Find top hub/auth scores per mode group (base R approach)
  top_hub <- df[df$hub.score == ave(df$hub.score, df$mode, FUN = function(x) max(x, na.rm = TRUE)), ]
  top_auth <- df[df$auth.score == ave(df$auth.score, df$mode, FUN = function(x) max(x, na.rm = TRUE)), ]

  # Plots
  p1 <- ggplot2::ggplot(df, ggplot2::aes(mode, wt.source, color = reorder(cells, cells))) +
    ggplot2::geom_point(size = df$hub.score * 10, alpha = 0.6) +
    ggplot2::coord_flip() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle("Outgoing Edgeweight") +
    ggplot2::geom_text(data = top_hub, ggplot2::aes(mode, wt.source, label = cells)) +
    ggplot2::ylab("Outgoing Edgeweight by Cell Type") +
    ggplot2::xlab("Signaling Family")
  if (normalize) {
    p1 <- p1 +
      ggplot2::ylab("Outgoing Edgeweight Fraction by Cell Type") +
      ggplot2::xlab("Signaling Family")
  }

  p2 <- ggplot2::ggplot(df, ggplot2::aes(mode, wt.sink, color = reorder(cells, cells))) +
    ggplot2::geom_point(size = df$auth.score * 10, alpha = 0.6) +
    ggplot2::coord_flip() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle("Incoming Edgeweight") +
    ggplot2::geom_text(data = top_auth, ggplot2::aes(mode, wt.sink, label = cells)) +
    ggplot2::ylab("Incoming Edgeweight by Cell Type") +
    ggplot2::xlab("Signaling Family")
  if (normalize) {
    p2 <- p2 +
      ggplot2::ylab("Incoming Edgeweight Fraction by Cell Type") +
      ggplot2::xlab("Signaling Family")
  }
  
  # Modify colors if desired
  if (!is.null(cols.use)) {
    p1 <- p1 + ggplot2::scale_colour_manual(values = cols.use)
    p2 <- p2 + ggplot2::scale_colour_manual(values = cols.use)
  }
  
  # Put legend on bottom
  legend <- cowplot::get_legend(
    p1 +
      ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 5))) +
      ggplot2::theme(legend.position = "bottom")
  )
  
  # Assemble plot
  plot.top <- cowplot::plot_grid(p1, p2, nrow = 1)
  return(cowplot::plot_grid(plot.top, legend, ncol = 1, rel_heights = c(1, .1)))
}

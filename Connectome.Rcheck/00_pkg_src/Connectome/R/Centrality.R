#' Centrality
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
#' @param group.by Default 'mode'. Determines how to subdivide network for centrality analysis.  Accepts 'mode','ligand','receptor', 'gene', or 'mechanism'.
#' @param ... Arguments passed to FilterConnectome
#' @return A ggplot object showing centrality analysis
#' @export

Centrality <- function(connectome,
                       cols.use = NULL,
                       weight.attribute = "weight_sc",
                       min.z = NULL,
                       normalize = TRUE,
                       group.by = "mode", ...) {
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
  ligands <- as.character(unique(master_sub$ligand))
  recepts <- as.character(unique(master_sub$receptor))
  genes <- as.character(union(unique(master_sub$ligand), unique(master_sub$receptor)))
  pairs <- as.character(unique(master_sub$pair))

  cells <- as.character(unique(union(master_sub$source, master_sub$target)))

  # Determine which type of network plot using switch for cleaner code
  groups <- switch(group.by,
    "mode" = modes,
    "ligand" = ligands,
    "receptor" = recepts,
    "gene" = genes,
    "mechanism" = pairs,
    modes  # default
  )

  # Process data
  n_groups <- length(groups)
  n_cells <- length(cells)
  result_list <- vector("list", n_groups * n_cells)
  idx <- 1

  for (i in seq_along(groups)) {
    # Subset data based on group type
    temp <- switch(group.by,
      "mode" = master_sub[master_sub$mode == groups[[i]], ],
      "ligand" = master_sub[master_sub$ligand == groups[[i]], ],
      "receptor" = master_sub[master_sub$receptor == groups[[i]], ],
      "gene" = master_sub[master_sub$ligand == groups[[i]] | master_sub$receptor == groups[[i]], ],
      "mechanism" = master_sub[master_sub$pair == groups[[i]], ],
      master_sub[master_sub$mode == groups[[i]], ]
    )

    # Skip if no data
    if (nrow(temp) == 0) next

    # Make network for analysis
    net <- igraph::graph_from_data_frame(temp, directed = TRUE)
    hub <- igraph::hub_score(net, weights = temp[, weight.attribute], scale = TRUE)$vector
    auth <- igraph::authority_score(net, weights = temp[, weight.attribute], scale = TRUE)$vector
    total.edgeweight <- sum(temp[, weight.attribute])

    # Get cumulative source/sink edgeweights per cell type
    for (j in seq_along(cells)) {
      wt.source <- sum(temp[temp$source == cells[[j]], weight.attribute])
      wt.sink <- sum(temp[temp$target == cells[[j]], weight.attribute])
      
      if (normalize && total.edgeweight > 0) {
        wt.source <- wt.source / total.edgeweight
        wt.sink <- wt.sink / total.edgeweight
      }
      
      # Store in pre-allocated list
      result_list[[idx]] <- data.frame(
        group = groups[[i]],
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

  # Combine all results
  df <- do.call(rbind, result_list[seq_len(idx - 1)])

  # Alphabetize
  df$group <- factor(df$group, levels = sort(as.character(unique(df$group)), decreasing = TRUE))

  # Plots
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = wt.source, color = as.factor(cells))) +
    ggplot2::geom_point(size = df$hub.score * 10, alpha = 0.6) +
    ggplot2::coord_flip() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle("Outgoing Centrality") +
    ggplot2::geom_text(data = df[df$hub.score == ave(df$hub.score, df$group, FUN = max), ], 
                       ggplot2::aes(group, wt.source, label = cells)) +
    ggplot2::ylab("Outgoing Edgeweight by Cell Type")
  if (normalize) {
    p1 <- p1 +
      ggplot2::ylab("Outgoing Edgeweight Fraction by Cell Type") +
      ggplot2::ylim(0, 1)
  }

  p2 <- ggplot2::ggplot(df, ggplot2::aes(group, wt.sink, color = as.factor(cells))) +
    ggplot2::geom_point(size = df$auth.score * 10, alpha = 0.6) +
    ggplot2::coord_flip() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle("Incoming Centrality") +
    ggplot2::geom_text(data = df[df$auth.score == ave(df$auth.score, df$group, FUN = max), ],
                       ggplot2::aes(group, wt.sink, label = cells)) +
    ggplot2::ylab("Incoming Edgeweight by Cell Type")
  if (normalize) {
    p2 <- p2 +
      ggplot2::ylab("Incoming Edgeweight Fraction by Cell Type") +
      ggplot2::ylim(0, 1)
  }

  # Define vertical axis label
  if (group.by == "mode") {
    p1 <- p1 + xlab("Network (by Family)")
    p2 <- p2 + xlab("Network (by Family)")
  }
  if (group.by == "ligand") {
    p1 <- p1 + xlab("Network (by Ligand)")
    p2 <- p2 + xlab("Network (by Ligand)")
  }
  if (group.by == "receptor") {
    p1 <- p1 + xlab("Network (by Receptor)")
    p2 <- p2 + xlab("Network (by Receptor)")
  }
  if (group.by == "gene") {
    p1 <- p1 + xlab("Network (by Gene)")
    p2 <- p2 + xlab("Network (by Gene)")
  }
  if (group.by == "mechanism") {
    p1 <- p1 + xlab("Network (by Mechanism)")
    p2 <- p2 + xlab("Network (by Mechanism)")
  }

  # Modify colors if desired
  if (!is.null(cols.use)) {
    p1 <- p1 + scale_colour_manual(values = cols.use)
    p2 <- p2 + scale_colour_manual(values = cols.use)
  }
  # Put legend on bottom
  legend <- get_legend(
    p1 +
      guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 5))) +
      theme(legend.position = "bottom")
  )
  # Assemble plot
  plot.top <- plot_grid(p1, p2, nrow = 1)
  return(plot_grid(plot.top, legend, ncol = 1, rel_heights = c(1, .1)))
}

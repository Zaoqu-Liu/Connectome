#' CompareCentrality
#'
#' Takes a list of connectomes and will compare sending- and receiving- centrality, side-by-side, for the given networks, or requested network subsets.
#'
#' @param connectome.list A named list of connectomic edgelists. The output plot will be split by the names of the list.
#' @param weight.attribute Column to use to define edgeweights for network analysis.  'weight_sc' or 'weight_norm'. Defaults to 'weight_norm'
#' @param min.z Minimum z-score for ligand and receptor.
#' @param ... Arguments passed to FilterConnectome.  Will be applied to each connectome within the list.
#' @param cols.use Desired colors for cell types, alphabetized. Defaults to standard ggplot colors.
#' @param normalize Default TRUE. Scales each mode to have equivalent x-axes.
#' @return A ggplot object comparing centrality across connectomes
#' @export

CompareCentrality <- function(connectome.list,
                              weight.attribute = "weight_norm",
                              min.z = NULL,
                              cols.use = NULL,
                              normalize = TRUE, ...) {
  # Check required packages
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required but not installed.")
  }
  if (!requireNamespace("ggthemes", quietly = TRUE)) {
    message("Package 'ggthemes' not found, using default theme.")
    use_ggthemes <- FALSE
  } else {
    use_ggthemes <- TRUE
  }

  if (weight.attribute == "weight_sc" && is.null(min.z)) {
    message("\nWeight attribute is 'weight_sc', recommend also setting min.z = 0 to avoid negative ligand and receptor scores")
  }

  # Add names to each connectome
  for (i in seq_along(connectome.list)) {
    connectome.list[[i]]$name <- names(connectome.list)[[i]]
  }

  # Filter all connectomes
  for (i in seq_along(connectome.list)) {
    connectome.list[[i]] <- FilterConnectome(connectome.list[[i]], min.z = min.z, remove.na = TRUE, ...)
  }

  # Process each connectome
  result_list <- vector("list", length(connectome.list))
  
  for (n in seq_along(connectome.list)) {
    master_func <- connectome.list[[n]]
    master_func$wt <- master_func[, weight.attribute]
    cells <- as.character(unique(union(master_func$source, master_func$target)))

    temp <- master_func
    net <- igraph::graph_from_data_frame(temp, directed = TRUE)
    hub <- igraph::hub_score(net, weights = temp$wt, scale = TRUE)$vector
    auth <- igraph::authority_score(net, weights = temp$wt, scale = TRUE)$vector
    total.edgeweight <- sum(temp[, weight.attribute])

    # Pre-allocate for inner loop
    cell_results <- vector("list", length(cells))
    
    for (j in seq_along(cells)) {
      wt.source <- sum(temp[temp$source == cells[[j]], "wt"])
      wt.sink <- sum(temp[temp$target == cells[[j]], "wt"])
      
      if (normalize && total.edgeweight > 0) {
        wt.source <- wt.source / total.edgeweight
        wt.sink <- wt.sink / total.edgeweight
      }
      
      cell_results[[j]] <- data.frame(
        cells = cells[[j]],
        hub.score = hub[cells[[j]]], 
        auth.score = auth[cells[[j]]],
        wt.source = wt.source, 
        wt.sink = wt.sink,
        name = master_func$name[1],
        stringsAsFactors = FALSE
      )
    }
    result_list[[n]] <- do.call(rbind, cell_results)
  }
  
  # Combine all results at once
  nn_total <- do.call(rbind, result_list)

  # Find top hub/auth scores per name group (base R approach)
  top_hub <- nn_total[nn_total$hub.score == ave(nn_total$hub.score, nn_total$name, FUN = max), ]
  top_auth <- nn_total[nn_total$auth.score == ave(nn_total$auth.score, nn_total$name, FUN = max), ]

  p1 <- ggplot2::ggplot(nn_total, ggplot2::aes(name, wt.source, color = as.factor(cells))) +
    ggplot2::geom_point(size = nn_total$hub.score * 10, alpha = 0.6) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 10))) +
    ggplot2::geom_text(data = top_hub, ggplot2::aes(name, wt.source, label = cells)) +
    ggplot2::ylim(0, max(nn_total$wt.source) * 1.1) +
    ggplot2::ggtitle("Outgoing Centrality") +
    ggplot2::ylab("Outgoing Edgeweight by Cell Type") +
    ggplot2::xlab("System") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)
    )
  
  if (use_ggthemes) {
    p1 <- p1 + ggthemes::theme_hc()
  }
  
  if (normalize) {
    p1 <- p1 + ggplot2::ylab("Outgoing Edgeweight Fraction by Cell Type")
  }
  
  p2 <- ggplot2::ggplot(nn_total, ggplot2::aes(name, wt.sink, color = as.factor(cells))) +
    ggplot2::geom_point(size = nn_total$auth.score * 10, alpha = 0.6) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 10))) +
    ggplot2::geom_text(data = top_auth, ggplot2::aes(name, wt.sink, label = cells)) +
    ggplot2::ylim(0, max(nn_total$wt.sink) * 1.1) +
    ggplot2::ggtitle("Incoming Centrality") +
    ggplot2::ylab("Incoming Edgeweight by Cell Type") +
    ggplot2::xlab("System") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)
    )
  
  if (use_ggthemes) {
    p2 <- p2 + ggthemes::theme_hc()
  }
  
  if (normalize) {
    p2 <- p2 + ggplot2::ylab("Incoming Edgeweight Fraction by Cell Type")
  }

  # Modify colors if desired
  if (!is.null(cols.use)) {
    p1 <- p1 + ggplot2::scale_colour_manual(values = cols.use)
    p2 <- p2 + ggplot2::scale_colour_manual(values = cols.use)
  }
  
  print(cowplot::plot_grid(p1, p2))
}

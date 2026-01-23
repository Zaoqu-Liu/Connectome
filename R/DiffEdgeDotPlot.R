#' DiffEdgeDotPlot
#'
#' Plotting function to make a DotPlot of differential edges for a network comparison. This function first finds all edges meeting the desired thresholding criteria, and then plots complete information regarding all mechanisms and celltype vectors implicated.
#'
#' @param differential.connectome A differential connectome object
#' @param sources.include Source nodes of interest. Output will be limited to edges coming from these sources.
#' @param targets.include Target nodes of interest. Output will be limited to edges landing on these targets.
#' @param features Gene of interest. Output will be limited to edges including these specific genes.
#' @param min.score Default NULL. Will limit output to edges with a differential score greater than this value.
#' @param min.pct Default NULL. Threshold to return clusterwise observations for both ligand and receptor. Only needs to be satisfied in connect.1 OR in connect.2.
#' @param verbose Whether to output feedback to user
#' @param infinity.to.max Default TRUE.  If TRUE, will create a pseudo value to replace values of "Inf"
#' @param ... Additional arguments (not currently used)
#' @return A ggplot object
#' @export

DiffEdgeDotPlot <- function(differential.connectome,
                            features = NULL,
                            sources.include = NULL,
                            targets.include = NULL,
                            min.score = NULL,
                            min.pct = NULL,
                            verbose = TRUE,
                            infinity.to.max = TRUE, ...) {
  # Check required packages - ggplot2 is in Imports so always available

  # Setup vector column
  differential.connectome$vector <- paste(differential.connectome$source,
    differential.connectome$target,
    sep = " - "
  )

  # Replace Inf/-Inf values with scaled max/min values
  if (infinity.to.max) {
    # Handle positive infinity
    pos_inf_lig <- is.infinite(differential.connectome$ligand.norm.lfc) & differential.connectome$ligand.norm.lfc > 0
    if (any(pos_inf_lig)) {
      differential.connectome$ligand.norm.lfc[pos_inf_lig] <- max(differential.connectome$ligand.norm.lfc[is.finite(differential.connectome$ligand.norm.lfc)], na.rm = TRUE) * 1.01
    }
    pos_inf_rec <- is.infinite(differential.connectome$recept.norm.lfc) & differential.connectome$recept.norm.lfc > 0
    if (any(pos_inf_rec)) {
      differential.connectome$recept.norm.lfc[pos_inf_rec] <- max(differential.connectome$recept.norm.lfc[is.finite(differential.connectome$recept.norm.lfc)], na.rm = TRUE) * 1.01
    }
    pos_inf_score <- is.infinite(differential.connectome$score) & differential.connectome$score > 0
    if (any(pos_inf_score)) {
      differential.connectome$score[pos_inf_score] <- max(differential.connectome$score[is.finite(differential.connectome$score)], na.rm = TRUE) * 1.01
    }

    # Handle negative infinity
    neg_inf_lig <- is.infinite(differential.connectome$ligand.norm.lfc) & differential.connectome$ligand.norm.lfc < 0
    if (any(neg_inf_lig)) {
      differential.connectome$ligand.norm.lfc[neg_inf_lig] <- min(differential.connectome$ligand.norm.lfc[is.finite(differential.connectome$ligand.norm.lfc)], na.rm = TRUE) * 1.01
    }
    neg_inf_rec <- is.infinite(differential.connectome$recept.norm.lfc) & differential.connectome$recept.norm.lfc < 0
    if (any(neg_inf_rec)) {
      differential.connectome$recept.norm.lfc[neg_inf_rec] <- min(differential.connectome$recept.norm.lfc[is.finite(differential.connectome$recept.norm.lfc)], na.rm = TRUE) * 1.01
    }
    neg_inf_score <- is.infinite(differential.connectome$score) & differential.connectome$score < 0
    if (any(neg_inf_score)) {
      differential.connectome$score[neg_inf_score] <- min(differential.connectome$score[is.finite(differential.connectome$score)], na.rm = TRUE) * 1.01
    }
  }

  # Filter
  data <- differential.connectome
  pre.filter <- nrow(data)

  # Subset based on min.score
  if (!is.null(min.score)) {
    data <- subset(data, score > min.score)
  }

  # Subset on nodes (cell types) of interest
  if (!is.null(sources.include)) {
    data <- subset(data, source %in% sources.include)
  }
  if (!is.null(targets.include)) {
    data <- subset(data, target %in% targets.include)
  }
  # Subset on features of interest
  if (!is.null(features)) {
    data <- subset(data, ligand %in% features | receptor %in% features)
  }
  # Subset based on min.pct
  if (!is.null(min.pct)) {
    data <- subset(data, pct.source.1 > min.pct | pct.source.2 > min.pct)
    data <- subset(data, pct.target.1 > min.pct | pct.target.2 > min.pct)
  }

  # Postfilter
  post.filter <- nrow(data)

  if (verbose) {
    message(paste("\nPre-filter edges: ", as.character(pre.filter)))
    message(paste("\nPost-filter edges: ", as.character(post.filter)))
    message("\nConnectome filtration completed")
  }


  # Identify vectors and mechanisms remaining
  vectors <- unique(data$vector)
  mechanisms <- unique(data$pair)

  # Subset to these
  data2 <- subset(differential.connectome, pair %in% mechanisms & vector %in% vectors)

  # Plot

  my_palette <- colorRampPalette(c("blue", "yellow", "red"), alpha = TRUE)(n = 399)


  p1 <- ggplot(data = data2, aes(
    x = vector, y = pair,
    # alpha = log(score+1)
  )) +
    geom_point(aes(size = log(score + 1), color = log(score + 1), stroke = 0)) +
    # theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_gradientn("log(score + 1)", colors = my_palette) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size = 14, colour = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.title = element_blank(),
      panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black")
    )

  return(p1)
}

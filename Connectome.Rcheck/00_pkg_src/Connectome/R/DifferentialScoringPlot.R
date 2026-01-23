#' DifferentialScoringPlot
#'
#' Currently in beta testing. Creates x3 aligned heatmaps allowing visualization of ligand, receptor, and perturbation scores for a given cell-system of interest.
#'
#' @param differential.connectome A differential connectome, made with DifferentialConnectome. May be filtered as desired prior to plotting.
#' @param sources.include Source nodes of interest. Output will be limited to edges coming from these sources.
#' @param targets.include Target nodes of interest. Output will be limited to edges landing on these targets.
#' @param features Gene of interest. Output will be limited to edges including these specific genes.
#' @param min.score Default NULL. Will limit output to edges with a differential score greater than this value.
#' @param min.pct Default NULL. Threshold to return clusterwise observations for both ligand and receptor. Only needs to be satisfied in connect.1 OR in connect.2.
#' @param verbose Whether to output feedback to user
#' @param infinity.to.max Default TRUE.  If TRUE, will create a pseudo value to replace values of "Inf"
#' @param aligned Default FALSE. If TRUE, will create edge-aligned heatmaps (duplicate rows and columns in first two plots, to dimension map all three plots)
#'
#' @export

DifferentialScoringPlot <- function(differential.connectome,
                                    features = NULL,
                                    sources.include = NULL,
                                    targets.include = NULL,
                                    min.score = NULL,
                                    min.pct = NULL,
                                    verbose = TRUE,
                                    infinity.to.max = TRUE,
                                    aligned = FALSE) {
  # ggplot2 and cowplot are in Imports

  # Setup vector column
  differential.connectome$vector <- paste(differential.connectome$source, differential.connectome$target, sep = " - ")

  # Save to data
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
  # Get total dataset for plotting (avoids grey squares)
  columns <- unique(data$vector)
  rows <- unique(data$pair)
  temp <- subset(differential.connectome, pair %in% rows & vector %in% columns)
  data <- temp

  # Replace Inf/-Inf values with scaled max/min values
  if (infinity.to.max) {
    # Handle positive infinity
    pos_inf_lig <- is.infinite(data$ligand.norm.lfc) & data$ligand.norm.lfc > 0
    if (any(pos_inf_lig)) {
      data$ligand.norm.lfc[pos_inf_lig] <- max(data$ligand.norm.lfc[is.finite(data$ligand.norm.lfc)], na.rm = TRUE) * 1.01
    }
    pos_inf_rec <- is.infinite(data$recept.norm.lfc) & data$recept.norm.lfc > 0
    if (any(pos_inf_rec)) {
      data$recept.norm.lfc[pos_inf_rec] <- max(data$recept.norm.lfc[is.finite(data$recept.norm.lfc)], na.rm = TRUE) * 1.01
    }
    pos_inf_score <- is.infinite(data$score) & data$score > 0
    if (any(pos_inf_score)) {
      data$score[pos_inf_score] <- max(data$score[is.finite(data$score)], na.rm = TRUE) * 1.01
    }

    # Handle negative infinity
    neg_inf_lig <- is.infinite(data$ligand.norm.lfc) & data$ligand.norm.lfc < 0
    if (any(neg_inf_lig)) {
      data$ligand.norm.lfc[neg_inf_lig] <- min(data$ligand.norm.lfc[is.finite(data$ligand.norm.lfc)], na.rm = TRUE) * 1.01
    }
    neg_inf_rec <- is.infinite(data$recept.norm.lfc) & data$recept.norm.lfc < 0
    if (any(neg_inf_rec)) {
      data$recept.norm.lfc[neg_inf_rec] <- min(data$recept.norm.lfc[is.finite(data$recept.norm.lfc)], na.rm = TRUE) * 1.01
    }
    neg_inf_score <- is.infinite(data$score) & data$score < 0
    if (any(neg_inf_score)) {
      data$score[neg_inf_score] <- min(data$score[is.finite(data$score)], na.rm = TRUE) * 1.01
    }
  }
  if (aligned) {
    # Set up replacement labels
    label.set <- unique(data[, c("ligand", "receptor", "pair")])
    label.set <- label.set[order(label.set$pair), ]
    label.set.2 <- unique(data[, c("source", "target", "vector")])
    label.set.2 <- label.set.2[order(label.set.2$vector), ]
    # Plot
    p1 <- ggplot(data, aes(x = vector, y = pair)) +
      geom_tile(aes(fill = ligand.norm.lfc)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
      ggtitle("Ligand Log2 Fold Change") +
      scale_y_discrete(labels = as.character(label.set$ligand)) +
      scale_x_discrete(labels = as.character(label.set.2$source)) +
      ylab("Ligand") +
      xlab("Source")

    p2 <- ggplot(data, aes(x = vector, y = pair)) +
      geom_tile(aes(fill = recept.norm.lfc)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
      ggtitle("Receptor Log2 Fold Change") +
      scale_y_discrete(labels = as.character(label.set$receptor)) +
      scale_x_discrete(labels = as.character(label.set.2$target)) +
      ylab("Receptor") +
      xlab("Target")

    p3 <- ggplot(data, aes(x = vector, y = pair)) +
      geom_tile(aes(fill = score)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
      ggtitle("Perturbation Score") +
      ylab("Mechanism") +
      xlab("Vector")

    plot_grid(p1, p2, p3, nrow = 1)
  } else {
    p1 <- ggplot(data, aes(x = source, y = ligand)) +
      geom_tile(aes(fill = ligand.norm.lfc)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
      ggtitle("Ligand Log2 Fold Change") +
      ylab("Ligand") +
      xlab("Source")

    p2 <- ggplot(data, aes(x = target, y = receptor)) +
      geom_tile(aes(fill = recept.norm.lfc)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
      ggtitle("Receptor Log2 Fold Change") +
      ylab("Receptor") +
      xlab("Target")

    p3 <- ggplot(data, aes(x = vector, y = pair)) +
      geom_tile(aes(fill = score)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
      ggtitle("Perturbation Score") +
      ylab("Mechanism") +
      xlab("Vector")

    plot_grid(p1, p2, p3, nrow = 1)
  }
}

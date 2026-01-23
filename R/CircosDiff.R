#' CircosDiff
#'
#' Currently in beta testing. A Circos Plot for a differential edgelist made using DifferentialConnectome. Note that negative log fold changes cannot be plotted; therefore, this plot uses the 'score' for the differential comparison, which is always positive in proportion to perturbation.
#'
#' @param differential.connectome A differential connectome, made with DifferentialConnectome. May be filtered as desired prior to plotting.
#' @param sources.include Source nodes of interest. Output will be limited to edges coming from these sources.
#' @param targets.include Target nodes of interest. Output will be limited to edges landing on these targets.
#' @param features Gene of interest. Output will be limited to edges including these specific genes.
#' @param min.score Default NULL. Will limit output to edges with a differential score greater than this value.
#' @param min.pct Default NULL. Threshold to return clusterwise observations for both ligand and receptor. Only needs to be satisfied in connect.1 OR in connect.2.
#' @param verbose Whether to output feedback to user
#' @param infinity.to.max Default TRUE.  If TRUE, will create a pseudo value to replace values of "Inf"
#' @param cols.use Optional. Colors for plotting nodes.
#' @param lab.cex Text size for gene names
#' @param edge.color.by.source Default TRUE - edges will be colored by their source cell type. If false, edges will be colored by receiving cell instead.
#' @param title Character string for title of plot.
#' @param small.gap Default 1. Amount of distance between sectors.  If the number of edges is very large, this will have to be reduced in size.
#' @param big.gap Default 10. Amount of distance between the source cells and the target cells (top and bottom arc of graph).  If the number of edges is very large, this can be reduced in size in addition to 'small.gap'


#' @export

CircosDiff <- function(differential.connectome,
                       features = NULL,
                       sources.include = NULL,
                       targets.include = NULL,
                       min.score = NULL,
                       min.pct = NULL,
                       verbose = TRUE,
                       infinity.to.max = TRUE,
                       edge.color.by.source = TRUE,
                       cols.use = NULL,
                       lab.cex = 1,
                       title = NULL,
                       small.gap = 1,
                       big.gap = 10) {
  # Check required packages
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required but not installed.")
  }
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required but not installed.")
  }

  data <- differential.connectome
  pre.filter <- nrow(data)

  # Setup vector column
  data$vector <- paste(data$source, data$target, sep = " - ")

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

  # CircosPlot Below:
  connectome <- data

  # Pull the dataframe of interest for plotting and format with weight as third column
  connectome$lig.stash <- as.character(connectome$ligand)
  connectome$rec.stash <- as.character(connectome$receptor)
  df <- data.frame(connectome %>% select(ligand, receptor))
  df$ligand <- make.unique(as.character(df$ligand))
  df$receptor <- make.unique(as.character(df$receptor))
  # df$weight <- connectome[,weight.attribute]
  temp <- connectome[, !colnames(connectome) %in% colnames(df)]
  df <- cbind(df, temp)

  # Squash ligands to single name if duplicates (different edges on same cell type)
  unique_ligs <- unique(df$lig.stash)
  for (i in seq_along(unique_ligs)) {
    lig_rows <- which(df$lig.stash == unique_ligs[i])
    temp <- df[lig_rows, ]
    unique_sources <- unique(temp$source)
    for (j in seq_along(unique_sources)) {
      source_rows <- lig_rows[temp$source == unique_sources[j]]
      dummy <- paste(rep(" ", j - 1), collapse = "")
      df$ligand[source_rows] <- paste0(as.character(df$lig.stash[source_rows[1]]), dummy)
    }
  }

  # Squash receptors back together to single name if they are duplicates (different edges on same cell type)
  unique_recs <- unique(df$rec.stash)
  for (i in seq_along(unique_recs)) {
    rec_rows <- which(df$rec.stash == unique_recs[i])
    temp <- df[rec_rows, ]
    unique_targets <- unique(temp$target)
    for (j in seq_along(unique_targets)) {
      target_rows <- rec_rows[temp$target == unique_targets[j]]
      dummy <- paste(rep(" ", j - 1), collapse = "")
      df$receptor[target_rows] <- paste0(as.character(df$rec.stash[target_rows[1]]), dummy)
    }
  }

  # Squash ligands back together, by cell type, if they are expressed on multiple cell types
  # temp <- subset(df,ligand != lig.stash) # this is a problem
  # if (nrow(temp)>0){
  #  for (i in 1:length(unique(temp$source))){
  #    temp2 <- subset(temp,source == unique(temp$source)[i])
  #    dummy <- paste(rep(' ',i),collapse = '') # Add number of spaces corresponding to number of unique sources
  #    df[rownames(temp2),]$ligand <- paste(as.character(temp2$lig.stash),dummy,sep='')
  #  }
  # }
  # Squash receptors back together, by cell type, if they are expressed on multiple cell types
  # temp <- subset(df,receptor != rec.stash) # this is a problem
  # if (nrow(temp)>0){
  #  for (i in 1:length(unique(temp$target))){
  #    temp2 <- subset(temp,target == unique(temp$target)[i])
  #    dummy <- paste(rep(' ',i),collapse = '') # Add number of spaces corresponding to number of unique targets
  #    df[rownames(temp2),]$receptor <- paste(as.character(temp2$rec.stash),dummy,sep='')
  #  }
  # }

  # Establish ordering (order) so that genes are grouped nicely by celltype
  source.order <- df[order(df$source), ]
  target.order <- df[order(df$target), ]
  source.order.un <- unique(source.order[, c("ligand", "source")])
  target.order.un <- unique(target.order[, c("receptor", "target")])

  source.order$id <- 1:nrow(source.order)
  target.order$id <- 1:nrow(target.order)
  source.order.un$id <- 1:nrow(source.order.un)
  target.order.un$id <- 1:nrow(target.order.un)

  sector.order.un <- c(
    as.character(source.order.un$ligand),
    as.character(target.order.un$receptor)
  )

  # Coloring setup
  if (is.null(cols.use)) {
    nodes <- as.character(unique(union(df$source, df$target)))
    cols.use <- hue_pal()(length(nodes))
    names(cols.use) <- nodes
    cols.use <- data.frame(cols.use)
    cols.use$cell <- rownames(cols.use)
  } else {
    cols.use <- data.frame(cols.use)
    cols.use$cell <- rownames(cols.use)
  }


  # Map to get ligand colorings (edges)
  map <- base::merge(source.order, cols.use, by.x = "source", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  lig.cols.edges <- as.character(map$cols.use)
  names(lig.cols.edges) <- map$ligand

  # Map to get receptor colorings (edges) # this does not work
  map <- base::merge(target.order, cols.use, by.x = "target", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  rec.cols.edges <- as.character(map$cols.use)
  names(rec.cols.edges) <- map$receptor

  # Map to get ligand colorings (sectors)
  map <- base::merge(source.order.un, cols.use, by.x = "source", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  lig.cols.sect <- as.character(map$cols.use)
  names(lig.cols.sect) <- map$ligand

  # Map to get receptor colorings (sectors)
  map <- base::merge(target.order.un, cols.use, by.x = "target", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  rec.cols.sect <- as.character(map$cols.use)
  names(rec.cols.sect) <- map$receptor

  # Make sector colors (grid.cols)
  sectors <- c(source.order.un$ligand, target.order.un$receptor)
  sector.cols <- c(as.character(lig.cols.sect), as.character(rec.cols.sect))

  # Plotting
  # Decide edge order and edge color order
  if (edge.color.by.source) {
    edge.color <- lig.cols.edges
    df.plot <- source.order
  } else {
    edge.color <- rec.cols.edges
    df.plot <- target.order
  }
  # Since differential, pull 'score' as the weight attribute (cannot plot negative values in log fold change data)

  df.plot <- df.plot[, c("ligand", "receptor", "score")]


  circos.clear()
  # circos.par(gap.degree = gap.degree)
  chordDiagram(df.plot,
    order = sector.order.un,
    col = edge.color,
    grid.col = sector.cols,
    directional = 1,
    direction.type = "arrows",
    link.arr.type = "big.arrow",
    annotationTrack = "grid",
    preAllocateTracks = 1,
    big.gap = big.gap,
    small.gap = small.gap
  )
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .01, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = lab.cex)
    # circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)
  # Make and add legend
  legend <- Legend(
    at = as.character(unique(union(df$source, df$target))),
    type = "grid",
    legend_gp = gpar(fill = as.character(cols.use[as.character(unique(union(df$source, df$target))), ]$cols.use)),
    title_position = "topleft",
    title = "Cell Type"
  )
  draw(legend, x = unit(20, "mm"), y = unit(20, "mm"), just = c("left", "bottom"))
  if (!is.null(title)) {
    title(title)
  }

  p1.base <- recordPlot()
  return(p1.base)
}

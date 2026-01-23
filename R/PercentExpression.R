#' Calculate Percent Expression
#'
#' Calculates the percentage of cells expressing each feature within each cluster.
#' This is an internal function used by CreateConnectome.
#'
#' @param object A Seurat object
#' @param assays Assay to use. Default "RNA"
#' @param features Features to calculate. Default NULL uses all features.
#' @param return.seurat Return a Seurat object instead of a list. Default FALSE.
#' @param add.ident Additional identity to include. Default NULL.
#' @param slot Data slot to use. Default "counts".
#' @param use.scale Deprecated. Use slot = "scale.data" instead.
#' @param use.counts Deprecated. Use slot = "counts" instead.
#' @param verbose Print progress messages. Default TRUE.
#' @param ... Additional arguments
#'
#' @return A list of data.frames with percent expression per cluster
#' @keywords internal

PercentExpression <- function(object,
                              assays = "RNA",
                              features = NULL,
                              return.seurat = FALSE,
                              add.ident = NULL,
                              slot = "counts",
                              use.scale = FALSE,
                              use.counts = FALSE,
                              verbose = TRUE, ...) {
  
  # Handle deprecated arguments

  if (use.scale) {
    .Deprecated(msg = "'use.scale' is deprecated, use slot = 'scale.data' instead")
    slot <- "scale.data"
  }
  if (use.counts) {
    .Deprecated(msg = "'use.counts' is deprecated, use slot = 'counts' instead")
    if (use.scale) {
      warning("Both 'use.scale' and 'use.counts' were set; using counts",
        call. = FALSE, immediate. = TRUE
      )
    }
    slot <- "counts"
  }

  # Get available assays
  object.assays <- Seurat::Assays(object = object)
  ident.orig <- Seurat::Idents(object = object)
  orig.levels <- levels(x = Seurat::Idents(object = object))
  ident.new <- c()

  if (!all(assays %in% object.assays)) {
    assays <- assays[assays %in% object.assays]
    if (length(assays) == 0) {
      stop("None of the requested assays are present in the object")
    } else {
      warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
    }
  }

  if (!is.null(x = add.ident)) {
    new.data <- Seurat::FetchData(object = object, vars = add.ident)
    new.ident <- paste(Seurat::Idents(object)[rownames(x = new.data)],
      new.data[, 1],
      sep = "_"
    )
    Seurat::Idents(object, cells = rownames(new.data)) <- new.ident
  }

  data.return <- list()

  for (i in seq_along(assays)) {
    if (verbose) {
      message("Processing assay: ", assays[i])
    }

    data.use <- Seurat::GetAssayData(object = object, assay = assays[i], slot = slot)
    features.assay <- features

    if (length(x = intersect(x = features, y = rownames(x = data.use))) < 1) {
      features.assay <- rownames(x = data.use)
    }

    # Get unique cluster identities
    clusters <- levels(x = Seurat::Idents(object))
    n_clusters <- length(clusters)

    # Ensure features are in data
    features.assay <- unique(x = intersect(x = features.assay, y = rownames(x = data.use)))
    n_features <- length(features.assay)

    if (verbose) {
      message("  Features: ", n_features, ", Clusters: ", n_clusters)
    }

    # Pre-allocate result matrix
    data.all <- matrix(NA_real_, nrow = n_features, ncol = n_clusters)
    rownames(data.all) <- features.assay
    colnames(data.all) <- clusters

    # Process each cluster
    for (j in seq_along(clusters)) {
      cluster_name <- clusters[j]
      temp.cells <- Seurat::WhichCells(object = object, idents = cluster_name)

      if (length(x = temp.cells) == 1) {
        data.temp <- as.numeric(data.use[features.assay, temp.cells] > 0)
        if (slot == "data") {
          data.temp <- as.numeric(expm1(data.use[features.assay, temp.cells]) > 0)
        }
      } else if (length(x = temp.cells) > 1) {
        data_subset <- data.use[features.assay, temp.cells, drop = FALSE]
        data.temp <- Matrix::rowSums(data_subset > 0) / length(temp.cells)
      } else {
        data.temp <- rep(0, n_features)
      }

      data.all[, j] <- data.temp

      if (verbose) {
        message("  Finished cluster ", cluster_name, " (", j, "/", n_clusters, ")")
      }

      if (i == 1) {
        if (length(temp.cells) > 0) {
          ident.new <- c(ident.new, as.character(x = ident.orig[temp.cells[1]]))
        }
      }
    }

    # Store results
    data.all <- as.data.frame(data.all)
    names(x = ident.new) <- levels(x = Seurat::Idents(object))
    data.return[[i]] <- data.all
    names(x = data.return)[i] <- assays[[i]]
  }

  if (return.seurat) {
    toRet <- Seurat::CreateSeuratObject(
      counts = data.return[[1]],
      project = "Average",
      assay = names(x = data.return)[1],
      ...
    )
    if (length(x = data.return) > 1) {
      for (i in seq(from = 2, to = length(data.return), by = 1)) {
        toRet[[names(x = data.return)[i]]] <- Seurat::CreateAssayObject(counts = data.return[[i]])
      }
    }
    if (Seurat::DefaultAssay(object = object) %in% names(x = data.return)) {
      Seurat::DefaultAssay(object = toRet) <- Seurat::DefaultAssay(object = object)
    }
    Seurat::Idents(toRet, cells = colnames(x = toRet)) <- ident.new[colnames(x = toRet)]
    Seurat::Idents(object = toRet) <- factor(
      x = Seurat::Idents(object = toRet),
      levels = as.character(x = orig.levels),
      ordered = TRUE
    )
    toRet <- Seurat::NormalizeData(object = toRet, verbose = verbose)
    toRet <- Seurat::ScaleData(object = toRet, verbose = verbose)
    return(toRet)
  } else {
    return(data.return)
  }
}

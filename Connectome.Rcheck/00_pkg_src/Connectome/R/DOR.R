#' DOR
#'
#' Calculate Diagnostic Odds Ratio for features. The DOR is a measure of gene
#' specificity for a cell cluster, defined as (TP/FP)/(FN/TN) where TP = cells
#' in cluster expressing gene, FP = cells outside cluster expressing gene,
#' FN = cells in cluster not expressing gene, TN = cells outside cluster not
#' expressing gene.
#'
#' @param object A Seurat object
#' @param ident Cell identity of interest
#' @param features Features of interest to calculate DOR
#' @param assay Assay to use for DOR calculation. Default "RNA"
#' @param pseudo Pseudocount to avoid division by zero. Default 0.5 (Haldane-Anscombe correction)
#' @return Named numeric vector of log-transformed DOR values
#' @keywords internal

DOR <- function(object,
                ident,
                features,
                assay = "RNA",
                pseudo = 0.5) {
  # Input validation
  if (!inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object")
  }
  
  cells <- Seurat::WhichCells(object = object, idents = ident)
  features <- as.character(features)

  # Get counts matrix using Seurat API (v4/v5 compatible)
  counts_matrix <- Seurat::GetAssayData(object, assay = assay, slot = "counts")

  # Ensure features are in matrix
  features <- intersect(features, rownames(counts_matrix))

  if (length(features) == 0) {
    warning("No valid features found in the counts matrix")
    return(numeric(0))
  }

  # Get all cells
  all_cells <- colnames(counts_matrix)
  cells_other <- setdiff(all_cells, cells)
  
  # Handle edge case where no cells match
  if (length(cells) == 0 || length(cells_other) == 0) {
    warning("No cells found for specified identity or no other cells available")
    result <- rep(NA_real_, length(features))
    names(result) <- features
    return(result)
  }

  # Extract relevant submatrices
  counts_in <- counts_matrix[features, cells, drop = FALSE]
  counts_out <- counts_matrix[features, cells_other, drop = FALSE]

  if (length(features) > 1) {
    # TP: cells in cluster expressing the gene
    TP <- Matrix::rowSums(counts_in > 0)

    # FN: cells in cluster NOT expressing the gene
    FN <- Matrix::rowSums(counts_in == 0)

    # FP: cells NOT in cluster expressing the gene
    FP <- Matrix::rowSums(counts_out > 0)

    # TN: cells NOT in cluster NOT expressing the gene
    TN <- Matrix::rowSums(counts_out == 0)
  } else {
    # Single feature
    TP <- sum(counts_in > 0)
    FN <- sum(counts_in == 0)
    FP <- sum(counts_out > 0)
    TN <- sum(counts_out == 0)
  }

  # Calculate DOR with Haldane-Anscombe correction (pseudo = 0.5)
  # This is the standard correction for 2x2 contingency tables
  # DOR = (TP / FP) / (FN / TN)
  DOR_val <- ((TP + pseudo) / (FP + pseudo)) / ((FN + pseudo) / (TN + pseudo))

  # Log transform
  DOR_val <- log(DOR_val)

  # Handle infinite values (can occur even with pseudocount in extreme cases)
  DOR_val[is.infinite(DOR_val)] <- NA_real_

  names(DOR_val) <- features
  return(DOR_val)
}

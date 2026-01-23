#' SplitMarkers
#'
#' Find differential markers between conditions within each cell type.
#' This function is specifically designed for comparing gene expression
#' between two conditions (e.g., organs) across cell types.
#'
#' @param object A Seurat object with cell type identities and a condition metadata column
#' @param GOI Character vector of genes of interest to test
#' @param condition.col Name of the metadata column containing condition information. Default "organ"
#' @param condition.1 First condition to compare. Default "liver"
#' @param condition.2 Second condition to compare. Default "lung"
#' @param ... Additional arguments passed to FindMarkers
#'
#' @return A data.frame with differential expression results including:
#'   gene, cell type, p-values, and comparison information
#' @keywords internal

SplitMarkers <- function(object, GOI, 
                         condition.col = "organ",
                         condition.1 = "liver",
                         condition.2 = "lung",
                         ...) {
  # Get cell types present in both conditions
  celltypes <- as.character(unique(Idents(object)))
  condition_table <- table(Idents(object), object@meta.data[[condition.col]])
  celltypes <- names(which(rowSums(condition_table != 0) == length(unique(object@meta.data[[condition.col]]))))

  if (length(celltypes) == 0) {
    warning("No cell types found in both conditions")
    return(data.frame())
  }

  celltypes.1 <- paste(celltypes, condition.1, sep = "_")
  celltypes.2 <- paste(celltypes, condition.2, sep = "_")

  object$celltype.condition <- paste(Idents(object), object@meta.data[[condition.col]], sep = "_")
  object$celltype <- Idents(object)
  Idents(object) <- "celltype.condition"

  # Filter genes to those present in object
  genes <- GOI[GOI %in% rownames(object)]
  
  if (length(genes) == 0) {
    warning("No genes of interest found in object")
    return(data.frame())
  }

  # Pre-allocate list for results
  result_list <- vector("list", length(celltypes))
  
  for (i in seq_along(celltypes)) {
    temp <- NULL
    try(
      temp <- FindMarkers(object,
        ident.1 = celltypes.1[i],
        ident.2 = celltypes.2[i],
        verbose = FALSE,
        features = genes,
        min.cells.group = 0, ...
      ),
      silent = TRUE
    )
    if (!is.null(temp)) {
      temp2 <- subset(temp, p_val_adj < 0.05)
      if (nrow(temp2) > 0) {
        temp2$gene <- rownames(temp2)
        temp2$cell <- celltypes[i]
        result_list[[i]] <- temp2
      }
    }
  }
  
  # Combine results
  diff.p <- do.call(rbind, result_list)
  
  if (is.null(diff.p) || nrow(diff.p) == 0) {
    return(data.frame())
  }
  
  diff.p$cell.gene <- paste(diff.p$cell, diff.p$gene, sep = ".")
  diff.p$Comparison <- paste(condition.1, "vs", condition.2, sep = ".")

  # Clean rownames

  rownames(diff.p) <- NULL
  return(diff.p)
}

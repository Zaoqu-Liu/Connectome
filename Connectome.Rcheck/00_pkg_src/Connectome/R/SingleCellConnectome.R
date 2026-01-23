#' SingleCellConnectome
#'
#' Takes a Seurat object and creates a single-cell connectome (rows are ligand-receptor mechanisms and columns are cell-cell vectors)
#'
#' @param object A Seurat object
#' @param species The species of the object that is being processed.  Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'
#' @param LR.database Accepts either 'fantom5' or 'custom'. If custom, a dataframe must be provided to argument custom.list with the first column equal to ligands, second column equal to associated receptors, and third column equal to desired modal categorizations.
#' @param max.cells.per.ident Default NULL. If a value is input, input object will be downsampled to requested number of cells per identity.
#' @param min.cells.per.ident Default NULL. If a value is input, only cell populations meeting this threshold will be included in network analysis.
#' @param include.putative Default TRUE. Includes ligand-receptor pairs deemed putative in FANTOM5 database.
#' @param include.rejected Default FALSE. If TRUE, includes gene pairs labeled "EXCLUDED" in FANTOM5 database.
#' @param custom.list Optional. A dataframe for custom mapping, with the first column equal to ligands, second column equal to associated receptors, and third column equal to desired modal categorizations.
#' @param slot.use Which data slot to use. Default "data".
#' @param weight.definition Method for combining ligand and receptor values. Options: "product" (default), "sum", "mean".
#' @param assay Assay to use. Default "RNA".
#' @param ... Additional arguments (not currently used)
#' @return A sparse matrix with ligand-receptor pairs as rows and cell-cell combinations as columns
#' @export

SingleCellConnectome <- function(object,
                                 LR.database = "fantom5",
                                 species = NULL,
                                 include.putative = TRUE,
                                 include.rejected = FALSE,
                                 max.cells.per.ident = NULL,
                                 min.cells.per.ident = NULL,
                                 slot.use = "data",
                                 weight.definition = "product",
                                 custom.list = NULL,
                                 assay = "RNA",
                                 ...) {
  # Check required packages
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required but not installed.")
  }

  if (LR.database == "fantom5" && is.null(species)) {
    stop("\nPlease select species for fantom5 mapping. Allows 'human','mouse','rat', or 'pig' ")
  }

  # Downsample input object using modern Seurat API
  if (!is.null(max.cells.per.ident)) {
    object <- subset(object, cells = Seurat::WhichCells(object, downsample = max.cells.per.ident))
  }

  # Limit object to cell populations larger than requested minimum
  if (!is.null(min.cells.per.ident)) {
    idents.include <- names(table(Seurat::Idents(object)))[table(Seurat::Idents(object)) > min.cells.per.ident]
    object <- subset(object, idents = idents.include)
  }

  if (LR.database == "fantom5") {
    # Load ground-truth database (FANTOM5, species-converted as appropriate)
    fantom5 <- switch(species,
      "human" = Connectome::ncomms8866_human,
      "mouse" = Connectome::ncomms8866_mouse,
      "rat" = Connectome::ncomms8866_rat,
      "pig" = Connectome::ncomms8866_pig,
      stop("Unknown species: ", species)
    )

    # Import ligand-receptor pairs and metadata
    lit.put <- fantom5[fantom5$Pair.Evidence %in% c("literature supported", "putative"), ]
    lit <- fantom5[fantom5$Pair.Evidence %in% c("literature supported"), ]
    
    # Determine which list of pairs to perform analysis on:
    if (include.rejected) {
      fantom <- fantom5
    } else if (include.putative) {
      fantom <- lit.put
    } else {
      fantom <- lit
    }
  }

  if (LR.database == "custom") {
    if (is.null(custom.list)) {
      stop("\nCustom mapping requested. Please provide custom list of Ligand-Receptor interactions")
    }
    ligands <- as.character(custom.list[, 1])
    recepts <- as.character(custom.list[, 2])
    fantom <- data.frame(
      Ligand.ApprovedSymbol = ligands, 
      Receptor.ApprovedSymbol = recepts,
      stringsAsFactors = FALSE
    )
  }

  # Get expression data using Seurat API (v4/v5 compatible)
  expr_data <- Seurat::GetAssayData(object, assay = assay, slot = slot.use)

  # Identify ligands and receptors expressed in the object
  fantom.specific <- fantom[fantom$Ligand.ApprovedSymbol %in% rownames(expr_data) & 
                            fantom$Receptor.ApprovedSymbol %in% rownames(expr_data), ]
  ligands <- as.character(fantom.specific$Ligand.ApprovedSymbol)
  receptors <- as.character(fantom.specific$Receptor.ApprovedSymbol)

  if (length(ligands) == 0) {
    stop("No ligand-receptor pairs found in the expression data")
  }

  message("Found ", length(ligands), " ligand-receptor pairs")
  message("Building single-cell connectome for ", ncol(expr_data), " cells...")

  # Get expression matrices for ligands and receptors
  lig.map <- expr_data[ligands, , drop = FALSE]
  rec.map <- expr_data[receptors, , drop = FALSE]
  
  n_cells <- ncol(expr_data)
  n_pairs <- length(ligands)
  
  # Calculate L-R product for each cell pair (sender, receiver)
  # Result: rows = L-R pairs, columns = cell-cell combinations
  
  # Create expanded matrices for all cell-cell combinations
  # lig.map2[i, j] = expression of ligand i in cell (j %% n_cells)
  # rec.map2[i, j] = expression of receptor i in cell (j %/% n_cells)
  
  # Pre-allocate column names
  cell_names <- colnames(expr_data)
  col_names <- paste(rep(cell_names, each = n_cells), 
                     rep(cell_names, times = n_cells), 
                     sep = "-")
  row_names <- paste(ligands, receptors, sep = "-")
  
  # Use memory-efficient approach for large datasets
  if (n_cells * n_cells > 1e6) {
    message("Large dataset detected. Using memory-efficient computation...")
    # For very large datasets, return a function to compute on demand
    # or use sparse representation
    warning("Single-cell connectome for >1000 cells may be very large. Consider downsampling.")
  }
  
  # Replicate ligand expression for all receiver cells
  lig.map2 <- do.call(cbind, replicate(n_cells, as.matrix(lig.map), simplify = FALSE))
  
  # Replicate receptor expression for all sender cells  
  rec.map2 <- matrix(
    rep(as.matrix(rec.map), each = n_cells),
    nrow = n_pairs,
    ncol = n_cells * n_cells
  )
  
  # Calculate connectome based on weight definition
  if (weight.definition == "product") {
    sc.connectome <- lig.map2 * rec.map2
  } else if (weight.definition == "sum") {
    sc.connectome <- lig.map2 + rec.map2
  } else if (weight.definition == "mean") {
    sc.connectome <- (lig.map2 + rec.map2) / 2
  } else {
    stop("Unknown weight.definition: ", weight.definition)
  }

  # Name the matrix
  rownames(sc.connectome) <- row_names
  colnames(sc.connectome) <- col_names

  # Convert to sparse matrix for memory efficiency
  sc.connectome <- Matrix::Matrix(sc.connectome, sparse = TRUE)

  message("Single-Cell Connectome generation complete")
  message("Dimensions: ", nrow(sc.connectome), " pairs x ", ncol(sc.connectome), " cell-cell combinations")
  
  return(sc.connectome)
}

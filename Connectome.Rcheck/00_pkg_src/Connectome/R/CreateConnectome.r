#' CreateConnectome
#'
#' Creates a connectomic edgelist from a Seurat object. The connectome represents
#' all possible ligand-receptor interactions between cell populations, weighted by
#' expression levels and statistical significance.
#'
#' @param object A Seurat object with normalized and scaled data
#' @param species Species for FANTOM5 database. Required if LR.database = 'fantom5'.
#'   Options: 'human', 'mouse', 'rat', or 'pig'
#' @param LR.database Ligand-receptor database to use. Either 'fantom5' (default) or
#'   'custom'. If 'custom', provide a dataframe via custom.list parameter.
#' @param max.cells.per.ident Maximum cells per identity for downsampling. Default NULL
#'   (no downsampling). Useful for large datasets.
#' @param min.cells.per.ident Minimum cells required per identity. Default NULL.
#'   Clusters below this threshold are excluded from analysis.
#' @param include.putative Include putative ligand-receptor pairs from FANTOM5.
#'   Default TRUE.
#' @param include.rejected Include pairs labeled "EXCLUDED" in FANTOM5. Default FALSE.
#' @param p.values Calculate Wilcoxon Rank Sum test p-values. Default TRUE.
#'   Set to FALSE for faster computation when statistical testing is not needed.
#' @param weight.definition.norm Method for edge weight from normalized data.
#'   Options: 'sum', 'mean', or 'product' (default).
#' @param weight.definition.scale Method for edge weight from scaled data.
#'   Options: 'sum', 'mean' (default), or 'product'.
#' @param custom.list Dataframe with custom ligand-receptor pairs. Required columns:
#'   ligand (col 1), receptor (col 2), mode/category (col 3).
#' @param calculate.DOR Calculate log-normalized Diagnostic Odds Ratio. Default FALSE.
#' @param assay Seurat assay to use. Default 'RNA'. Also accepts 'SCT'.
#' @param parallel Use parallel processing for p-value calculations. Default FALSE.
#' @param n.cores Number of cores for parallel processing. Default: detectCores() - 1
#' @param ... Additional arguments passed to internal functions
#'
#' @return A data.frame containing the connectomic edgelist with columns:
#'   source, target, ligand, receptor, pair, mode, expression values,
#'   z-scores, percent expressing, edge weights, and optionally p-values and DOR.
#'
#' @details
#' The function calculates edge weights using the formula specified by
#' weight.definition.norm and weight.definition.scale parameters:
#' \itemize{
#'   \item 'product': ligand * receptor (captures co-expression)
#'   \item 'sum': ligand + receptor (additive contribution)
#'   \item 'mean': (ligand + receptor) / 2 (average contribution)
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' connectome <- CreateConnectome(seurat_obj, species = "human")
#'
#' # Without p-value calculation (faster)
#' connectome <- CreateConnectome(seurat_obj, species = "human", p.values = FALSE)
#'
#' # With custom ligand-receptor list
#' my_lr <- data.frame(ligand = c("IL6", "TNF"), receptor = c("IL6R", "TNFRSF1A"),
#'                     mode = c("cytokine", "cytokine"))
#' connectome <- CreateConnectome(seurat_obj, LR.database = "custom", custom.list = my_lr)
#' }
#'
#' @export

CreateConnectome <- function(object,
                             LR.database = "fantom5",
                             species = NULL,
                             include.putative = TRUE,
                             include.rejected = FALSE,
                             p.values = TRUE,
                             max.cells.per.ident = NULL,
                             min.cells.per.ident = NULL,
                             weight.definition.norm = "product",
                             weight.definition.scale = "mean",
                             custom.list = NULL,
                             calculate.DOR = FALSE,
                             assay = "RNA",
                             parallel = FALSE,
                             n.cores = NULL, ...) {
  # Check required packages (R CMD check compliant)
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  if (LR.database == "fantom5" && is.null(species)) {
    stop("\nPlease select species for fantom5 mapping. Allows 'human','mouse','rat', or 'pig' ")
  }

  # Downsample input object using modern Seurat API
  if (!is.null(max.cells.per.ident)) {
    message("Downsampling to ", max.cells.per.ident, " cells per identity...")
    object <- subset(object, cells = Seurat::WhichCells(object, downsample = max.cells.per.ident))
  }

  # Limit object to cell populations larger than requested minimum
  if (!is.null(min.cells.per.ident)) {
    message("Filtering cell populations with >", min.cells.per.ident, " cells...")
    idents.include <- names(table(Seurat::Idents(object)))[table(Seurat::Idents(object)) > min.cells.per.ident]
    object <- subset(object, idents = idents.include)
  }

  # Load ligand-receptor database
  if (LR.database == "fantom5") {
    message("Loading FANTOM5 database for species: ", species)

    if (species == "human") {
      fantom5 <- Connectome::ncomms8866_human
    } else if (species == "mouse") {
      fantom5 <- Connectome::ncomms8866_mouse
    } else if (species == "rat") {
      fantom5 <- Connectome::ncomms8866_rat
    } else if (species == "pig") {
      fantom5 <- Connectome::ncomms8866_pig
    }

    # Import ligand-receptor pairs and metadata
    lit.put <- fantom5[fantom5$Pair.Evidence %in% c("literature supported", "putative"), ]
    lit <- fantom5[fantom5$Pair.Evidence %in% c("literature supported"), ]

    # Determine which list of pairs to perform analysis on:
    if (include.putative) {
      ligands <- as.character(lit.put[, 2])
      recepts <- as.character(lit.put[, 4])
      modes <- as.character(lit.put[, "mode"])
    } else {
      ligands <- as.character(lit[, 2])
      recepts <- as.character(lit[, 4])
      modes <- as.character(lit[, "mode"])
    }
    if (include.rejected) {
      ligands <- as.character(fantom5[, 2])
      recepts <- as.character(fantom5[, 4])
      modes <- as.character(fantom5[, "mode"])
    }
  }

  if (LR.database == "custom") {
    if (is.null(custom.list)) {
      stop("\nCustom mapping requested. Please provide custom list of Ligand-Receptor interactions")
    }
    ligands <- as.character(custom.list[, 1])
    recepts <- as.character(custom.list[, 2])
    modes <- as.character(custom.list[, 3])
  }

  message("Total L-R pairs in database: ", length(ligands))

  # Identify ligands and receptors expressed in the object (Seurat v4/v5 compatible)
  assay_data <- Seurat::GetAssayData(object, assay = assay, slot = "data")
  ligands.use <- intersect(ligands, rownames(assay_data))
  recepts.use <- intersect(recepts, rownames(assay_data))
  genes.use <- union(ligands.use, recepts.use)

  message(
    "Expressed genes found: ", length(genes.use),
    " (", length(ligands.use), " ligands, ", length(recepts.use), " receptors)"
  )

  # Create averages and other relevant cluster-wise metrics
  message("Calculating cluster-wise expression metrics...")
  cluster.avgs <- Seurat::AverageExpression(object, features = genes.use, assays = assay)[[assay]]
  cluster.avgs.scale <- Seurat::AverageExpression(object, features = genes.use, slot = "scale.data", assays = assay)[[assay]]

  cluster.pcts <- PercentExpression(object, features = genes.use, slot = "counts", assays = assay, verbose = FALSE)[[assay]]

  if (calculate.DOR) {
    message("Calculating Diagnostic Odds Ratios...")
    nodes <- as.character(names(table(Seurat::Idents(object))))
    cluster.DORs <- data.frame(row.names = genes.use)
    for (i in seq_along(nodes)) {
      temp <- DOR(object, ident = nodes[i], features = genes.use)
      temp <- data.frame(temp)
      colnames(temp) <- nodes[i]
      cluster.DORs <- cbind(cluster.DORs, temp)
    }
  }

  # Convert to data.table for efficient merging
  cluster.avgs.dt <- data.table::as.data.table(cluster.avgs, keep.rownames = "gene")
  cluster.avgs.scale.dt <- data.table::as.data.table(cluster.avgs.scale, keep.rownames = "gene")
  cluster.pcts.dt <- data.table::as.data.table(cluster.pcts, keep.rownames = "gene")

  if (calculate.DOR) {
    cluster.DORs.dt <- data.table::as.data.table(cluster.DORs, keep.rownames = "gene")
  }

  # Create ligand and receptor lookup tables
  lr.pairs.dt <- data.table::data.table(
    id = seq_along(ligands),
    ligand = ligands,
    receptor = recepts,
    mode = modes
  )

  # Merge with expression data
  ligands.dt <- merge(lr.pairs.dt[, c("id", "ligand", "receptor", "mode")],
    cluster.avgs.dt,
    by.x = "ligand", by.y = "gene",
    all.x = TRUE, sort = FALSE
  )
  ligands.dt <- merge(ligands.dt,
    cluster.avgs.scale.dt,
    by.x = "ligand", by.y = "gene",
    all.x = TRUE, sort = FALSE,
    suffixes = c("", ".scale")
  )
  ligands.dt <- merge(ligands.dt,
    cluster.pcts.dt,
    by.x = "ligand", by.y = "gene",
    all.x = TRUE, sort = FALSE,
    suffixes = c("", ".pct")
  )

  recepts.dt <- merge(lr.pairs.dt[, c("id", "receptor")],
    cluster.avgs.dt,
    by.x = "receptor", by.y = "gene",
    all.x = TRUE, sort = FALSE
  )
  recepts.dt <- merge(recepts.dt,
    cluster.avgs.scale.dt,
    by.x = "receptor", by.y = "gene",
    all.x = TRUE, sort = FALSE,
    suffixes = c("", ".scale")
  )
  recepts.dt <- merge(recepts.dt,
    cluster.pcts.dt,
    by.x = "receptor", by.y = "gene",
    all.x = TRUE, sort = FALSE,
    suffixes = c("", ".pct")
  )

  if (calculate.DOR) {
    ligands.dt <- merge(ligands.dt,
      cluster.DORs.dt,
      by.x = "ligand", by.y = "gene",
      all.x = TRUE, sort = FALSE,
      suffixes = c("", ".DOR")
    )
    recepts.dt <- merge(recepts.dt,
      cluster.DORs.dt,
      by.x = "receptor", by.y = "gene",
      all.x = TRUE, sort = FALSE,
      suffixes = c("", ".DOR")
    )
  }

  # Sort by id to maintain original order
  data.table::setorder(ligands.dt, id)
  data.table::setorder(recepts.dt, id)

  # Calculate p-values if required
  if (p.values) {
    message("Calculating p-values using Wilcoxon Rank Sum test...")

    if (parallel && requireNamespace("future", quietly = TRUE) &&
      requireNamespace("future.apply", quietly = TRUE)) {
      message("Using parallel processing for p-value calculations...")

      if (is.null(n.cores)) {
        n.cores <- max(1, parallel::detectCores() - 1)
      }

      old_plan <- future::plan()
      future::plan(future::multisession, workers = n.cores)

      cluster.p.values <- Seurat::FindAllMarkers(object,
        assay = assay,
        features = genes.use,
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = 0,
        return.thresh = Inf
      )

      future::plan(old_plan)
    } else {
      cluster.p.values <- Seurat::FindAllMarkers(object,
        assay = assay,
        features = genes.use,
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = 0,
        return.thresh = Inf
      )
    }

    cluster.p.values$source.ligand <- paste(cluster.p.values$cluster, cluster.p.values$gene, sep = " - ")
    cluster.p.values$receptor.target <- paste(cluster.p.values$gene, cluster.p.values$cluster, sep = " - ")
    cluster.p.values.dt <- data.table::as.data.table(cluster.p.values)
  }

  # Generate full connectome
  message("Generating connectome...")
  sources <- colnames(cluster.avgs)
  targets <- colnames(cluster.avgs)
  n_sources <- length(sources)
  n_targets <- length(targets)
  n_pairs <- nrow(lr.pairs.dt)
  n_total <- n_sources * n_targets * n_pairs

  message(
    "Creating ", n_sources, " x ", n_targets, " = ", n_sources * n_targets,
    " cell-cell connections with ", n_pairs, " L-R pairs each"
  )
  message("Total edges: ", format(n_total, big.mark = ","))

  # Pre-allocate list for edge data
  connectome_list <- vector("list", n_sources * n_targets)
  idx <- 1

  # Progress bar
  pb <- txtProgressBar(min = 0, max = n_sources, style = 3)

  for (i in seq_len(n_sources)) {
    source_name <- sources[i]

    for (j in seq_len(n_targets)) {
      target_name <- targets[j]

      # Extract relevant columns for this source-target pair
      lig_scale_cols <- grep(paste0("^", source_name, "\\.scale$"), names(ligands.dt), value = TRUE)
      lig_pct_cols <- grep(paste0("^", source_name, "\\.pct$"), names(ligands.dt), value = TRUE)

      rec_scale_cols <- grep(paste0("^", target_name, "\\.scale$"), names(recepts.dt), value = TRUE)
      rec_pct_cols <- grep(paste0("^", target_name, "\\.pct$"), names(recepts.dt), value = TRUE)

      # Build edge data
      edge_data <- data.table::data.table(
        source = source_name,
        target = target_name,
        ligand = ligands.dt$ligand,
        receptor = ligands.dt$receptor,
        pair = paste(ligands.dt$ligand, ligands.dt$receptor, sep = " - "),
        mode = ligands.dt$mode,
        ligand.expression = ligands.dt[[source_name]],
        recept.expression = recepts.dt[[target_name]],
        ligand.scale = if (length(lig_scale_cols) > 0) ligands.dt[[lig_scale_cols[1]]] else ligands.dt[[paste0(source_name, ".scale")]],
        recept.scale = if (length(rec_scale_cols) > 0) recepts.dt[[rec_scale_cols[1]]] else recepts.dt[[paste0(target_name, ".scale")]],
        percent.source = if (length(lig_pct_cols) > 0) ligands.dt[[lig_pct_cols[1]]] else ligands.dt[[paste0(source_name, ".pct")]],
        percent.target = if (length(rec_pct_cols) > 0) recepts.dt[[rec_pct_cols[1]]] else recepts.dt[[paste0(target_name, ".pct")]]
      )

      if (calculate.DOR) {
        lig_dor_cols <- grep(paste0("^", source_name, "\\.DOR$"), names(ligands.dt), value = TRUE)
        rec_dor_cols <- grep(paste0("^", target_name, "\\.DOR$"), names(recepts.dt), value = TRUE)

        edge_data$DOR.source <- if (length(lig_dor_cols) > 0) ligands.dt[[lig_dor_cols[1]]] else ligands.dt[[paste0(source_name, ".DOR")]]
        edge_data$DOR.target <- if (length(rec_dor_cols) > 0) recepts.dt[[rec_dor_cols[1]]] else recepts.dt[[paste0(target_name, ".DOR")]]
      }

      connectome_list[[idx]] <- edge_data
      idx <- idx + 1
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Combine all edge data
  message("\nCombining edge data...")
  connectome <- data.table::rbindlist(connectome_list, use.names = TRUE, fill = TRUE)

  # Log transform the expression values (since they were averaged outside of log space)
  message("Log-transforming expression values...")
  connectome$ligand.expression <- log1p(connectome$ligand.expression)
  connectome$recept.expression <- log1p(connectome$recept.expression)

  # Add weights (normalized slot)
  message("Calculating edge weights...")
  if (weight.definition.norm == "sum") {
    connectome$weight_norm <- connectome$ligand.expression + connectome$recept.expression
  } else if (weight.definition.norm == "mean") {
    connectome$weight_norm <- (connectome$ligand.expression + connectome$recept.expression) / 2
  } else if (weight.definition.norm == "product") {
    connectome$weight_norm <- connectome$ligand.expression * connectome$recept.expression
  } else {
    message("\nNo appropriate parameter specified for weight.definition.norm")
  }

  # Add weights (scaled slot)
  if (weight.definition.scale == "sum") {
    connectome$weight_sc <- connectome$ligand.scale + connectome$recept.scale
  } else if (weight.definition.scale == "mean") {
    connectome$weight_sc <- (connectome$ligand.scale + connectome$recept.scale) / 2
  } else if (weight.definition.scale == "product") {
    connectome$weight_sc <- connectome$ligand.scale * connectome$recept.scale
  } else {
    message("\nNo appropriate parameter specified for weight.definition.scale")
  }

  # Additional columns
  message("Adding metadata columns...")
  connectome$vector <- paste(connectome$source, connectome$target, sep = " - ")
  connectome$edge <- paste(connectome$source, connectome$ligand, connectome$receptor, connectome$target, sep = " - ")
  connectome$source.ligand <- paste(connectome$source, connectome$ligand, sep = " - ")
  connectome$receptor.target <- paste(connectome$receptor, connectome$target, sep = " - ")
  if (!is.null(species)) {
    connectome$species <- species
  }

  # Add p-values if required
  if (p.values) {
    message("Mapping p-values to edges...")
    
    # Prepare p-value lookup tables
    pval_lig <- cluster.p.values.dt[, c("source.ligand", "p_val_adj")]
    names(pval_lig)[2] <- "p_val_adj.lig"
    
    pval_rec <- cluster.p.values.dt[, c("receptor.target", "p_val_adj")]
    names(pval_rec)[2] <- "p_val_adj.rec"

    # Merge ligand p-values
    connectome <- merge(connectome, pval_lig,
      by = "source.ligand",
      all.x = TRUE,
      sort = FALSE
    )

    # Merge receptor p-values
    connectome <- merge(connectome, pval_rec,
      by = "receptor.target",
      all.x = TRUE,
      sort = FALSE
    )

    # Remove rows with non-sensical mappings
    connectome <- connectome[!is.na(connectome$source) & !is.na(connectome$target), ]
  }

  # Reorganize columns for presentation
  first_cols <- c("source", "target", "ligand", "receptor", "pair", "mode")
  other_cols <- setdiff(names(connectome), first_cols)
  data.table::setcolorder(connectome, c(first_cols, other_cols))

  # Convert back to data.frame for compatibility with original function
  connectome <- as.data.frame(connectome)

  message(paste("\nConnectome generation complete"))

  return(connectome)
}

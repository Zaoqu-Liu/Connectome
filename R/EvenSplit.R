#' EvenSplit
#'
#' Splits a Seurat object using SplitObject, but makes sure that for each identity, an even number of cells are sampled.
#'
#' @param object A Seurat object
#' @param split.by Parameter by which to split the object (e.g., a metadata column name)
#' @return A named list of Seurat objects, one for each level of split.by
#' @export

EvenSplit <- function(object, split.by) {
  # Number of cells grabb-able from each condition, per cell type
  nums <- table(Idents(object), object@meta.data[[split.by]])
  mins <- apply(nums, 1, FUN = min)
  split.object <- SplitObject(object, split.by = split.by)

  # Equalize the grabs across condition
  split.down <- list()
  for (i in seq_along(split.object)) {
    pile <- list()
    for (j in seq_along(mins)) {
      temp <- subset(split.object[[i]], idents = names(mins[j]))
      temp <- subset(temp, cells = WhichCells(temp, downsample = mins[j]))
      pile[[j]] <- temp
    }
    split.down[[i]] <- pile
  }
  # Squash into a list of objects
  for (i in seq_along(split.down)) {
    if (length(split.down[[i]]) > 1) {
      split.down[[i]] <- merge(split.down[[i]][[1]], split.down[[i]][2:length(split.down[[i]])])
    } else {
      split.down[[i]] <- split.down[[i]][[1]]
    }
  }
  names(split.down) <- names(split.object)
  return(split.down)
}

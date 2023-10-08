#' @title Conversion to AtgnatData
#'
#' @concept atgnat-object
#'
#' @rdname as.AtgnatData
#' @param x An object to take counts from
#' @param ... additional arguments passed to object-specific methods.
#'
#' @return An [AtgnatData] object
#'
#' @export
setGeneric(name = "as.AtgnatData",
           signature = c("x"),
           def = function(x, ...) standardGeneric("as.AtgnatData"))

#' @rdname as.AtgnatData
#'
#' @import methods
#' @import Matrix
#'
#' @export
setMethod(
  f = "as.AtgnatData",
  signature = c(x="data.frame"),
  definition = function(x){
    return(methods::new("AtgnatData", pseudobulks = x, bins = colnames(x)))
  }
)


#' @rdname as.AtgnatData
#' @param pseudotime_slot The [SingleCellExperiment::SingleCellExperiment] slot containing the pseudotime values for each cell
#' @param n_bins The number of bins to create.
#' @param split_by The technique used to split the bins. The default `pseudotime_range` picks the
#' bin for a cell based on a constant range of pseudotime. `cells` picks the bin for a cell based on
#' an even number of cells per bin.
#'
#' @import methods
#'
#' @export
setMethod(
  f = "as.AtgnatData",
  signature = c(x="SingleCellExperiment"),
  definition = function(x, pseudotime_slot="slingPseudotime_1", n_bins=20, split_by="pseudotime_range"){

    pseudotime_sce = assign_pseduotime_bins(x, split_by="pseudotime_range", n_bins=n_bins, pseudotime_slot=pseudotime_slot)
    c = SingleCellExperiment::normcounts(pseudotime_sce)

    bin_ids = sort(unique(pseudotime_sce$pseudotime_bin))
    pseudobulks = as.data.frame(Matrix::rowSums(SingleCellExperiment::normcounts(pseudotime_sce)))
    for (i in bin_ids) {
      pb = as.data.frame(Matrix::rowSums(SingleCellExperiment::normcounts(subset(pseudotime_sce, , pseudotime_sce@colData[["pseudotime_bin"]] == i))))
      colnames(pb) = i
      pseudobulks = cbind(pseudobulks, pb)
    }
    pseudobulks = pseudobulks[2:length(pseudobulks)]

    return(methods::new("AtgnatData", pseudobulks = pseudobulks, bins = bin_ids))
  }
)


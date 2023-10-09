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
#' @param pseudotime_slot The [SingleCellExperiment::SingleCellExperiment] slot
#' containing pseudotime values for each cell to be passed to [assign_pseudotime_bins()].
#' @param n_bins The number of bins to create, passed to [assign_pseudotime_bins()].
#' @param split_by The split_by method to be passed on to [assign_pseudotime_bins()].
#'
#' @import methods
#'
#' @export
setMethod(
  f = "as.AtgnatData",
  signature = c(x="SingleCellExperiment"),
  definition = function(x, pseudotime_slot="slingPseudotime_1", n_bins=20, split_by="pseudotime_range"){

    pseudotime_sce = assign_pseudotime_bins(x, split_by, n_bins=n_bins, pseudotime_slot=pseudotime_slot)
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


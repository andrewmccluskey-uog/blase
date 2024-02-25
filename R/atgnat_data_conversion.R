#' @title Conversion to AtgnatData
#'
#' @concept blase-object
#'
#' @rdname as.AtgnatData
#' @param x An object to take counts from
#' @param ... additional arguments passed to object-specific methods.
#'
#' @return An [AtgnatData] object
#'
#' @export
#' @inherit AtgnatData-class examples
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
#'
#' @import methods
#' @import Matrix
#'
#' @import rlang
#'
#' @export
setMethod(
  f = "as.AtgnatData",
  signature = c(x="Seurat"),
  definition = function(x, pseudotime_slot="slingPseudotime_1", n_bins=20, split_by="pseudotime_range"){
    rlang::check_installed("Seurat", reason = "to handle Seurat objects.")
    sce = Seurat::as.SingleCellExperiment(x)
    return(as.AtgnatData(sce, pseudotime_slot=pseudotime_slot, n_bins=n_bins, split_by=split_by))
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

    # bin_ids = sort(unique(pseudotime_sce$pseudotime_bin)) historically
    bin_ids = sort(unique(pseudotime_sce@colData[["pseudotime_bin"]]))
    pseudobulks = list()

    for (i in bin_ids) {
      bin_subset_sce = subset(pseudotime_sce, , pseudotime_bin == i)
      # TODO below line does not work but we generate a note with the alternative (above).
      # pseudotime_sym = ggplot2::sym("pseudotime_bin")
      # bin_subset_sce = subset(pseudotime_sce, , {pseudotime_sym} == i)
      counts = SingleCellExperiment::normcounts(bin_subset_sce)
      pseudobulks[[i]] = SingleCellExperiment::normcounts(bin_subset_sce)
    }
    return(methods::new("AtgnatData", pseudobulk_bins = pseudobulks, bins = bin_ids))
  }
)


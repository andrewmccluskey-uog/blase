#' @title Conversion to AtgnatData
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

    if (is.na(match(split_by, c("pseudotime_range", "cells")))) {
      stop("split_by must be 'pseudotime_range' or 'cells'")
    }

    # TODO add check for pseudotime_slot existing

    pseudotime_sce = subset(x, , !is.na(x@colData[pseudotime_slot]))
    pseudotime = pseudotime_sce@colData[[pseudotime_slot]]

    c = SingleCellExperiment::normcounts(pseudotime_sce)

    if (split_by == "pseudotime_range") {
      min_pdt = 0
      max_pdt = ceiling(max(pseudotime))

      bin_size = max_pdt/n_bins
      bin_upper_limits = seq(bin_size, max_pdt, by=bin_size)

      ## Put cells into the right bins
      pseudotime_sce$pseudotime_bin = ceiling(pseudotime / bin_size)

      pseudotime_sce$pseudotime_bin[pseudotime_sce$pseudotime_bin==0] = 1
    } else {
      pseudotime_order = order(pseudotime, decreasing = FALSE)

      ncells = ncol(SingleCellExperiment::normcounts(pseudotime_sce))
      cells_per_bin = floor(ncells/n_bins)
      pseudotime_ordered_cells = rownames(pseudotime_sce@colData)[pseudotime_order]

      pseudotime_sce$pseudotime_bin = n_bins

      for (i in seq_len(n_bins)) {
        cells_for_bin = pseudotime_ordered_cells[(i*cells_per_bin-cells_per_bin+1) : (i*cells_per_bin)]
        pseudotime_sce[,cells_for_bin]$pseudotime_bin = i
      }
    }

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

#' @title Assign Pseudotime Bins to a source object's metadata
#'
#' @concept atgnat-object
#'
#' @rdname assign_pseduotime_bins
#' @param x An object to add metadata to
#' @param ... additional arguments passed to object-specific methods.
#'
#' @export
setGeneric(name = "assign_pseduotime_bins",
           signature = c(x="x"),
           def = function(x, split_by="pseudotime_range", n_bins=20, ...) standardGeneric("assign_pseduotime_bins"))

#' @rdname assign_pseduotime_bins
#'
#' @export
setMethod(
  f = "assign_pseduotime_bins",
  signature = c(x="SingleCellExperiment"),
  definition = function(x, split_by, n_bins, pseudotime_slot="slingPseudotime_1"){

    if (is.na(match(split_by, c("pseudotime_range", "cells")))) {
      stop("split_by must be 'pseudotime_range' or 'cells'")
    }

    pseudotime_sce = subset(x, , !is.na(x@colData[pseudotime_slot]))
    pseudotime = pseudotime_sce@colData[[pseudotime_slot]]

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

    return(pseudotime_sce)


  }
)

#' @rdname assign_pseduotime_bins
#'
#' @export
setMethod(
  f = "assign_pseduotime_bins",
  signature = c(x="data.frame"),
  definition = function(x, split_by, n_bins){
    stop("Can't update bulk data, using each sample as bins.")
  }
)

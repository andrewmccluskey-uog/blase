#' Create pseudotime bins
#'
#' @description
#' Creates the pseudotime bins that will be used to characterise pseudobulked pseudotime.
#'
#'
#' @param pseudotime_sce The SingleCellExperiment object to add the bins metadata to
#' @param number_of_bins The number of bins to be created. Defaults to 20
#' @param pseudotime_slot The slot to use for pseudotime values. Defaults to "slingPseudotime_1"
#' @param split_by The technique used to split the bins. The default "pseudotime_range" picks the
#' bin for a cell based on a constant range of pseudotime. "cells" picks the bin for a cell based on
#' an even number of cells per bin.
#'
#' @return The SCE object with the bins metadata created in the "pseudotime_bins" slot.
#' @export
#'
#'
#' @examples
#' library(SingleCellExperiment, quiet=TRUE)
#' library(atgnat)
#' counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
#' sce <- SingleCellExperiment::SingleCellExperiment(counts)
#' sce$pseudotime = seq_len(10)
#' sce = create_pseudotime_bins(sce, 5, pseudotime_slot="pseudotime")
#' sce$pseudotime_bins

## -----------------------------------------------------------------------------
create_pseudotime_bins <- function(pseudotime_sce, number_of_bins=20, pseudotime_slot="slingPseudotime_1", split_by="pseudotime_range") {

  if (is.na(match(split_by, c("pseudotime_range", "cells")))) {
    stop("split_by must be 'pseudotime_range' or 'cells'")
  }

  ## Remove NA pseudotimes
  pseudotime_sce = subset(pseudotime_sce, , !is.na(pseudotime_sce@colData[pseudotime_slot]))

  pseudotime = pseudotime_sce@colData[[pseudotime_slot]]

  if (split_by == "cells") {
    pseudotime_sce = PRIVATE_apply_pseudotime_even_cells(pseudotime_sce, number_of_bins, pseudotime)
  } else {
    pseudotime_sce = PRIVATE_apply_pseudotime_even_pseudotime_range(pseudotime_sce, number_of_bins, pseudotime)
  }

  return(pseudotime_sce)
}

PRIVATE_apply_pseudotime_even_pseudotime_range <- function(pseudotime_sce, number_of_bins, pseudotime) {

  min_pdt = 0
  max_pdt = ceiling(max(pseudotime[!is.na(pseudotime)]))

  bin_size = max_pdt/number_of_bins
  bin_upper_limits = seq(bin_size, max_pdt, by=bin_size)

  ## Put cells into the right bins
  pseudotime_sce$pseudotime_bin = ceiling(pseudotime / bin_size)

  pseudotime_sce$pseudotime_bin[pseudotime_sce$pseudotime_bin==0] = 1

  return(pseudotime_sce)

}


PRIVATE_apply_pseudotime_even_cells = function(pseudotime_sce, number_of_bins, pseudotime) {
  pseudotime_order = order(pseudotime, decreasing = FALSE)

  ncells = ncol(SingleCellExperiment::counts(pseudotime_sce))
  cells_per_bin = floor(ncells/number_of_bins)
  pseudotime_ordered_cells = rownames(pseudotime_sce@colData)[pseudotime_order]

  pseudotime_sce$pseudotime_bin = number_of_bins

  for (i in seq_len(number_of_bins)) {
    cells_for_bin = pseudotime_ordered_cells[(i*cells_per_bin-cells_per_bin+1) : (i*cells_per_bin)]
    pseudotime_sce[,cells_for_bin]$pseudotime_bin = i
  }

  return(pseudotime_sce)
}

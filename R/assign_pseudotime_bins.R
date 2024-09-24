#' @title Assign Pseudotime Bins to a source object's metadata
#'
#' @concept blase-object
#'
#' @rdname assign_pseudotime_bins
#' @param x An object to add metadata to.
#' @param split_by The technique used to split the bins. The default
#' `pseudotime_range` picks the bin for a cell based on a constant range of
#' pseudotime. `cells` picks the bin for a cell based on an even number of
#' cells per bin.
#' @param n_bins The number of bins to split the cells into.
#' @param pseudotime_slot The [SingleCellExperiment::SingleCellExperiment]
#' slot containing the pseudotime values for each cell.
#' @param ... For arguments passed to other functions. Unused.
#'
#' @returns A copy of x where cells are annotated with their
#' pseudotime bin.
#'
#' @export
#' @inherit MappingResult-class examples
setGeneric(
    name = "assign_pseudotime_bins",
    signature = c(x = "x"),
    def = function(x,
                   split_by = "pseudotime_range",
                   n_bins = 20,
                   pseudotime_slot = "slingPseudotime_1",
                   ...) {
        standardGeneric("assign_pseudotime_bins")
    }
)

#' @rdname assign_pseudotime_bins
#'
#' @export
setMethod(
    f = "assign_pseudotime_bins",
    signature = c(x = "SingleCellExperiment"),
    definition = function(x, split_by, n_bins,
                          pseudotime_slot = "slingPseudotime_1") {
        PRIVATE_assign_pseudotime_bins_validate_inputs(
            x, split_by,
            pseudotime_slot
        )

        pseudotime_sce <- subset(x, , !is.na(
            SummarizedExperiment::colData(x)[pseudotime_slot]
        ))
        pseudotime <- (
            SummarizedExperiment::colData(pseudotime_sce)[[pseudotime_slot]]
        )

        if (split_by == "pseudotime_range") {
            min_pdt <- 0
            max_pdt <- ceiling(max(pseudotime))

            bin_size <- max_pdt / n_bins
            bin_upper_limits <- seq(bin_size, max_pdt, by = bin_size)

            ## Put cells into the right bins
            pseudotime_sce$pseudotime_bin <- ceiling(pseudotime / bin_size)

            pseudotime_sce$pseudotime_bin[
                pseudotime_sce$pseudotime_bin == 0
            ] <- 1
        } else {
            pseudotime_order <- order(pseudotime, decreasing = FALSE)

            ncells <- ncol(SingleCellExperiment::normcounts(pseudotime_sce))
            cells_per_bin <- floor(ncells / n_bins)
            pseudotime_ordered_cells <- rownames(pseudotime_sce@colData)[
                pseudotime_order
            ]

            pseudotime_sce$pseudotime_bin <- n_bins

            for (i in seq_len(n_bins)) {
                cells_for_bin <- pseudotime_ordered_cells[
                    (i * cells_per_bin - cells_per_bin + 1):(i * cells_per_bin)
                ]
                pseudotime_sce[, cells_for_bin]$pseudotime_bin <- i
            }
        }

        return(pseudotime_sce)
    }
)

#' @rdname assign_pseudotime_bins
#'
#' @export
setMethod(
    f = "assign_pseudotime_bins",
    signature = c(x = "data.frame"),
    definition = function(x, split_by, n_bins, pseudotime_slot = "slingPseudotime_1") {
        stop("Can't update bulk data, using each sample as bins.")
    }
)

#' @rdname assign_pseudotime_bins
#'
#' @import rlang
#'
#' @export
setMethod(
    f = "assign_pseudotime_bins",
    signature = c(x = "Seurat"),
    definition = function(x,
                          split_by,
                          n_bins,
                          pseudotime_slot = "slingPseudotime_1") {
        rlang::check_installed("Seurat", reason = "to handle Seurat objects.")
        sce <- Seurat::as.SingleCellExperiment(x)
        sce <- assign_pseudotime_bins(sce, split_by, n_bins, pseudotime_slot)
        x$pseudotime_bin <- sce$pseudotime_bin
        return(x)
    }
)

PRIVATE_assign_pseudotime_bins_validate_inputs <- function(
    x,
    split_by,
    pseudotime_slot) {
    if (is.na(match(split_by, c("pseudotime_range", "cells")))) {
        stop("split_by must be 'pseudotime_range' or 'cells'")
    }

    if (!any(colnames(x@colData) == pseudotime_slot)) {
        stop("Pseudotime slot '", pseudotime_slot, "' does not exist")
    }

    # TODO check for rownames and colnames existing
}

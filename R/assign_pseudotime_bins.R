#' @title Assign Pseudotime Bins to a source object's metadata
#'
#' @concept blase-object
#'
#' @rdname assign_pseudotime_bins
#' @param x An object to add metadata to.
#' @param split_by String. The technique used to split the bins. The default
#' `pseudotime_range` picks the bin for a cell based on a constant range of
#' pseudotime. `cells` picks the bin for a cell based on an even number of
#' cells per bin.
#' @param n_bins Integer. The number of bins to split the cells into.
#' @param pseudotime_slot String or Vector of Strings. The name of the
#' [SingleCellExperiment::SingleCellExperiment] slot containing the
#' pseudotime values for each cell.
#' @param ... For arguments passed to other functions. Unused.
#'
#' @returns A copy of x where cells are annotated with their
#' pseudotime bin.
#'
#' @export
#' @inherit MappingResult examples
setGeneric(
    name = "assign_pseudotime_bins",
    signature = c(x = "x"),
    def = function(x, split_by = "pseudotime_range", n_bins = 20,
                   pseudotime_slot = "slingPseudotime_1", ...) {
        standardGeneric("assign_pseudotime_bins")
    }
)

#' @rdname assign_pseudotime_bins
#'
#' @importFrom SummarizedExperiment colData
#'
#' @export
setMethod(
    f = "assign_pseudotime_bins",
    signature = c(x = "SingleCellExperiment"),
    definition = function(x, split_by, n_bins,
                          pseudotime_slot = "slingPseudotime_1") {
      PRIVATE_assign_pseudotime_bins_validate_inputs(
          x, split_by, n_bins, pseudotime_slot)

      # order lineages
      if (length(pseudotime_slot) > 1) {
        # Order lineages by length (largest to smallest max value)
        get_lineage_max = function(lineage) {
          maxval <- max(SummarizedExperiment::colData(x)[[lineage]], na.rm = TRUE)
          return(maxval)
        }
        lineage_lengths = vapply(pseudotime_slot, get_lineage_max, double(1), USE.NAMES = FALSE)

        lineages_long_to_short = pseudotime_slot[sort.list(lineage_lengths, decreasing=TRUE)]
        # Assign cells to each lineage
        x$original_lineage <- NA
        x$combined_pseudotime <- -1
        for (lin in lineages_long_to_short) {
          pseudotime <- SummarizedExperiment::colData(x)[[lin]]
          cell_mask = is.na(x$original_lineage) & !is.na(pseudotime)

          cell_order = rank(pseudotime[cell_mask])
          previous_max = max(0, max(x$combined_pseudotime))

          x$combined_pseudotime[cell_mask] <- previous_max + cell_order

          x$original_lineage[cell_mask] <- lin
        }

        # remove cells with no lineage

        # Based on cell counts per lineage, assign number of bins
        lineage_cell_counts = table(x$original_lineage)

        nominal_bin_size = floor(ncol(x[,x$combined_pseudotime!=-1]) / n_bins)
        message(paste("Nominal cells per bin:", nominal_bin_size))

        get_lineage_bins = function(lineage) {
          n_cells = lineage_cell_counts[[lineage]]
          if (n_cells < nominal_bin_size){
            return(1)
          }
          return(floor(n_cells/nominal_bin_size))
        }
        lineage_n_bins = vapply(rev(lineages_long_to_short), get_lineage_bins, numeric(1))

        x$pseudotime_bin = NA
        for (lin in lineages_long_to_short) {
          bins_for_lineage = lineage_n_bins[[lin]]

          y = x[, x$original_lineage==lin & !is.na(x$original_lineage)]
          pdt = SummarizedExperiment::colData(y)[[lin]]

          SummarizedExperiment::colData(y)[[lin]] = pdt - min(pdt, na.rm=T)
          y = PRIVATE_assign_pseudotime_bins_for_one_lineage(y, split_by, bins_for_lineage, lin)

          # include 0 as an option as for first sample, all are NA, so max()
          # returns -Inf
          x$pseudotime_bin = replace(
            x$pseudotime_bin,
            x$original_lineage==lin & !is.na(x$original_lineage),
            y$pseudotime_bin + max(x$pseudotime_bin, 0, na.rm = TRUE))
        }

        return(x)

      } else {
        return(PRIVATE_assign_pseudotime_bins_for_one_lineage(x, split_by, n_bins, pseudotime_slot))
      }

    }
)

backup = function() {
  # I think we may want to do this slightly differently so that we can
  # enforce a hard cutoff between lineage end/starts.

  # Check if more than one slot given.
  if (length(pseudotime_slot) > 1) {
    # Order lineages by length (largest to smallest max value)
    get_linage_max = function(lineage) {
      maxval <- max(SummarizedExperiment::colData(x)[[lineage]], na.rm = TRUE)
      return(maxval)
    }
    lineage_lengths = vapply(pseudotime_slot, get_linage_max, double(1), USE.NAMES = FALSE)
    x$combined_bin_pseudotime <- -1
    x$original_lineage <- "?"

    # for each lineage
    #   Cells without a combined_bin_pseudotime value are given one
    #   of their order.

    for (lin in pseudotime_slot[sort.list(lineage_lengths, decreasing=TRUE)]) {
      pseudotime <- (
        SummarizedExperiment::colData(x)[[lin]]
      )
      cell_mask = x$combined_bin_pseudotime==-1 & !is.na(pseudotime)
      cell_order = rank(pseudotime[cell_mask])
      previous_max = max(0, max(x$combined_bin_pseudotime))

      x$combined_bin_pseudotime[cell_mask] <- previous_max + cell_order
      x$original_lineage[cell_mask] <- lin
    }
    x$combined_bin_pseudotime = x$combined_bin_pseudotime - 1

    pseudotime_slot = "combined_bin_pseudotime"
  }
  return(PRIVATE_assign_pseudotime_bins_for_one_lineage(x, split_by, n_bins, pseudotime_slot))
}

#' @rdname assign_pseudotime_bins
#'
#' @export
setMethod(
    f = "assign_pseudotime_bins",
    signature = c(x = "data.frame"),
    definition = function(x, split_by, n_bins,
                          pseudotime_slot = "slingPseudotime_1") {
        stop("Can't update bulk data, using each sample as bins.")
    }
)

#' @rdname assign_pseudotime_bins
#'
#' @importFrom rlang check_installed
#' @importFrom Seurat as.SingleCellExperiment
#'
#' @export
setMethod(
    f = "assign_pseudotime_bins",
    signature = c(x = "Seurat"),
    definition = function(x, split_by, n_bins,
                          pseudotime_slot = "slingPseudotime_1") {
        rlang::check_installed("Seurat", reason = "to handle Seurat objects.")
        sce <- Seurat::as.SingleCellExperiment(x)
        sce <- assign_pseudotime_bins(sce, split_by, n_bins, pseudotime_slot)
        x$pseudotime_bin <- sce$pseudotime_bin
        return(x)
    }
)

#' @keywords internal
#' @importFrom SummarizedExperiment colData
PRIVATE_assign_pseudotime_bins_validate_inputs <- function(
    x,
    split_by,
    n_bins,
    pseudotime_slot) {
    if (is.na(match(split_by, c("pseudotime_range", "cells")))) {
        stop("split_by must be 'pseudotime_range' or 'cells'")
    }

    if (n_bins<=0) {
      stop("Number of bins must be greater than 0")
    }

    if (n_bins > ncol(x)) {
      stop("Number of bins must be less than cells in object")
    }

    if (length(pseudotime_slot) == 1) {
      if (!any(colnames(SummarizedExperiment::colData(x)) == pseudotime_slot)) {
        stop("Pseudotime slot '", pseudotime_slot, "' does not exist")
      }
    } else {
      if (split_by=="pseudotime_range") {
        stop("Cannot use pseudotime range assignment for ",
             "multiple trajectories.")
      }
    }

    for (s in pseudotime_slot) {
      if (!any(colnames(SummarizedExperiment::colData(x)) == s)) {
        stop("Pseudotime slot '", s, "' does not exist")
      }
    }

}

#' @keywords internal
PRIVATE_assign_pseudotime_bins_for_one_lineage <- function(
    x,
    split_by,
    n_bins,
    pseudotime_slot) {

  pseudotime_sce <- subset(x, , !is.na(
    SummarizedExperiment::colData(x)[pseudotime_slot]
  ))
  pseudotime <- (
    SummarizedExperiment::colData(pseudotime_sce)[[pseudotime_slot]]
  )
  if (min(pseudotime) != 0) {
    stop("pseudotime must start at 0")
  }

  if (split_by == "pseudotime_range") {
    min_pdt <- 0
    max_pdt <- max(pseudotime)

    bin_size <- (max_pdt) / n_bins
    bin_upper_limits <- seq(bin_size, max_pdt, by = bin_size)

    ## Put cells into the right bins
    pseudotime_sce$pseudotime_bin <- ceiling(pseudotime / bin_size)

    ## In slingshot there is one cell which starts exactly at 0 which
    ## which we want to include in bin 1. We can't access a bin in R
    ## which has the name 0 as R is 1 indexed.
    pseudotime_sce$pseudotime_bin[
      pseudotime_sce$pseudotime_bin == 0
    ] <- 1
  } else {
    ncells <- ncol(SingleCellExperiment::normcounts(pseudotime_sce))
    cells_per_bin <- floor(ncells / n_bins)
    pseudotime_ordered_cells <- rownames(
      SummarizedExperiment::colData(pseudotime_sce))[
        order(pseudotime, decreasing = FALSE)
      ]

    pseudotime_sce$pseudotime_bin <- n_bins

    for (i in seq_len(n_bins)) {
      cells_for_bin <- pseudotime_ordered_cells[
        (i * cells_per_bin - cells_per_bin + 1):(i * cells_per_bin)
      ]
      # pseudotime_sce[, cells_for_bin]$pseudotime_bin <- i

      print(i)
      pseudotime_sce$pseudotime_bin = replace(pseudotime_sce$pseudotime_bin, colnames(pseudotime_sce) %in% cells_for_bin, i)
    }
  }

  return(pseudotime_sce)

}

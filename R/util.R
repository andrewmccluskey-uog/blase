#' Get a pseudobulk of bins with at least 2 replicates
#'
#' This function will try to create a pseudobulked count matrix for the bins.
#' When a replicate has too few cells, it is discounted. If only one exists,
#' then we sample from it twice to create the pseudobulks.
#'
#' @concept util
#'
#' @param pseudotime_sce The Single Cell Experiment object to
#' get the bins from
#' @param min_cells_for_bulk The minimum cells to look for per
#' replicate and bin.
#' @param replicate_slot The slot in the Single Cell Experiment
#' that contains replicate information
#'
#' @return A dataframe pseudobulk counts matrix.
#' @export
#'
#' @examples
#' library(SingleCellExperiment, quietly = TRUE)
#' library(blase)
#' counts <- matrix(rpois(1000, lambda = 10), ncol = 100, nrow = 10)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(normcounts = counts, counts = counts / 2)
#' )
#' sce$pseudotime <- seq_len(100)
#' colnames(sce) <- seq_len(100)
#' rownames(sce) <- as.character(seq_len(10))
#' sce <- assign_pseudotime_bins(sce,
#'     n_bins = 5,
#'     pseudotime_slot = "pseudotime", split_by = "cells"
#' )
#' sce$replicate <- rep(c(1, 2), 50)
#' result <- get_bins_as_bulk(
#'     sce,
#'     min_cells_for_bulk = 1,
#'     replicate_slot = "replicate"
#' )
#' result
get_bins_as_bulk <- function(
    pseudotime_sce, min_cells_for_bulk = 50,
    replicate_slot = "replicate") {
    output <- data.frame()
    for (bin_id in seq_len(max(pseudotime_sce$pseudotime_bin))) {
        bin_specific_sce <- subset(
            pseudotime_sce, ,
            pseudotime_sce@colData[["pseudotime_bin"]] == bin_id
        )

        counts <- table(bin_specific_sce@colData[[replicate_slot]])

        replicates_with_more_than_minimum <- rownames(
            as.data.frame(counts[counts > min_cells_for_bulk])
        )

        if (length(replicates_with_more_than_minimum) >= 2) {
            pseudobulks <- PRIVATE_get_bins_as_bulk_bin_has_replicates(
                replicates_with_more_than_minimum,
                bin_specific_sce,
                replicate_slot,
                bin_id
            )
        } else if (length(replicates_with_more_than_minimum) == 1) {
            pseudobulks <- PRIVATE_get_bins_as_bulk_bin_without_replicates(
                replicates_with_more_than_minimum,
                bin_specific_sce,
                replicate_slot,
                bin_id
            )
        } else {
            message(
                "Couldn't create pseudobulks due to",
                " too few cells in every replicate for bin ",
                bin_id
            )
            next()
        }

        if (ncol(output) == 0) {
            output <- pseudobulks
        } else {
            output <- merge(output, pseudobulks, by = "row.names")
            rownames(output) <- output$Row.names
            output <- output[, c(-1)]
        }
    }

    return(output)
}

PRIVATE_get_bins_as_bulk_bin_has_replicates <- function(
    replicates_with_more_than_minimum,
    bin_specific_sce,
    replicate_slot,
    bin_id) {
    pseudobulks <- data.frame()
    for (rep_id in replicates_with_more_than_minimum) {
        bin_specific_rep_specific_sce_pseudobulk <- as.data.frame(
            Matrix::rowSums(SingleCellExperiment::counts(subset(
                bin_specific_sce, ,
                bin_specific_sce@colData[[replicate_slot]] == rep_id
            )))
        )

        if (ncol(pseudobulks) == 0) {
            pseudobulks <- bin_specific_rep_specific_sce_pseudobulk
            colnames(pseudobulks) <- c(
                paste0("bin_", bin_id, "_rep_", rep_id)
            )
        } else {
            pseudobulks <- merge(
                x = pseudobulks,
                y = bin_specific_rep_specific_sce_pseudobulk,
                by = "row.names"
            )
            rownames(pseudobulks) <- pseudobulks$Row.names
            pseudobulks <- pseudobulks[, c(-1)]
            colnames(pseudobulks) <- append(
                colnames(pseudobulks)[seq_len(ncol(pseudobulks) - 1)],
                paste0("bin_", bin_id, "_rep_", rep_id)
            )
        }
    }
    return(pseudobulks)
}

PRIVATE_get_bins_as_bulk_bin_without_replicates <- function(
    replicates_with_more_than_minimum,
    bin_specific_sce,
    replicate_slot,
    bin_id) {
    bin_specific_pseudobulk <- subset(
        bin_specific_sce, ,
        bin_specific_sce@colData[[replicate_slot]] ==
            replicates_with_more_than_minimum[1]
    )

    cell_count_to_sample <- ceiling((ncol(
        SingleCellExperiment::counts(bin_specific_pseudobulk)
    ) * 0.75))
    counts_for_bulk <- SingleCellExperiment::counts(
        bin_specific_pseudobulk
    )
    pseudobulks <- merge(
        x = as.data.frame(rowSums(counts_for_bulk[
            , sample(ncol(counts_for_bulk), size = cell_count_to_sample)
        ])),
        y = as.data.frame(rowSums(counts_for_bulk[
            , sample(ncol(counts_for_bulk), size = cell_count_to_sample)
        ])),
        by = "row.names"
    )
    rownames(pseudobulks) <- pseudobulks$Row.names
    pseudobulks <- pseudobulks[, c(-1)]
    colnames(pseudobulks) <- c(
        paste0(
            "bin_",
            bin_id,
            "_rep_",
            replicates_with_more_than_minimum[1],
            "_1"
        ),
        paste0(
            "bin_",
            bin_id,
            "_rep_",
            replicates_with_more_than_minimum[1],
            "_2"
        )
    )
    return(pseudobulks)
}

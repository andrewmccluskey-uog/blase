#' @title Conversion to BlaseData
#'
#' @concept blase-object
#'
#' @rdname as.BlaseData
#' @param x An object to take counts from
#' @param ... additional arguments passed to object-specific methods.
#'
#' @return An [BlaseData] object
#'
#' @export
#' @inherit BlaseData-class examples
setGeneric(
    name = "as.BlaseData",
    signature = c("x"),
    def = function(x, ...) standardGeneric("as.BlaseData")
)

#' @rdname as.BlaseData
#'
#' @import methods
#' @import Matrix
#'
#' @export
setMethod(
    f = "as.BlaseData",
    signature = c(x = "data.frame"),
    definition = function(x) {
        return(methods::new("BlaseData", pseudobulks = x, bins = colnames(x)))
    }
)

#' @rdname as.BlaseData
#'
#' @import methods
#' @import Matrix
#'
#' @import rlang
#'
#' @export
setMethod(
    f = "as.BlaseData",
    signature = c(x = "Seurat"),
    definition = function(x, pseudotime_slot = "slingPseudotime_1",
                          n_bins = 20, split_by = "pseudotime_range") {
        rlang::check_installed("Seurat", reason = "to handle Seurat objects.")
        sce <- Seurat::as.SingleCellExperiment(x)
        return(
            as.BlaseData(
                sce,
                pseudotime_slot = pseudotime_slot,
                n_bins = n_bins,
                split_by = split_by
            )
        )
    }
)


#' @rdname as.BlaseData
#' @param pseudotime_slot The [SingleCellExperiment::SingleCellExperiment]
#' slot containing pseudotime values for each cell to be passed to
#' [assign_pseudotime_bins()].
#' @param n_bins The number of bins to create, passed to
#' [assign_pseudotime_bins()].
#' @param split_by The split_by method to be passed on to
#' [assign_pseudotime_bins()].
#'
#' @import methods
#'
#' @export
setMethod(
    f = "as.BlaseData",
    signature = c(x = "SingleCellExperiment"),
    definition = function(x, pseudotime_slot = "slingPseudotime_1",
                          n_bins = 20, split_by = "pseudotime_range") {
        pseudotime_sce <- assign_pseudotime_bins(
            x,
            split_by,
            n_bins = n_bins,
            pseudotime_slot = pseudotime_slot
        )

        bin_ids <- sort(unique(pseudotime_sce@colData[["pseudotime_bin"]]))
        pseudobulks <- list()

        for (i in bin_ids) {
            bin_subset_sce <- pseudotime_sce[
                , SingleCellExperiment::colData(pseudotime_sce)[[
                    "pseudotime_bin"
                ]] == i
            ]
            counts <- SingleCellExperiment::normcounts(bin_subset_sce)
            pseudobulks[[i]] <- SingleCellExperiment::normcounts(bin_subset_sce)
        }
        return(methods::new(
            "BlaseData",
            pseudobulk_bins = pseudobulks,
            bins = bin_ids
        ))
    }
)

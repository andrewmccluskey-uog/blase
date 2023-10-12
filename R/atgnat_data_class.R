#' Atgnat Data Object
#'
#' For creation details, see [as.AtgnatData()]
#'
#' @concept atgnat-object
#'
#' @slot pseudobulks data.frame. Each column is a timepoint sample and each row is a gene.
#' @slot bins list. A list of bin names for each timepoint.
#' @slot genes list. A list of the genes selected for discriminating timepoints.
#'
#' @return An [AtgnatData] object
#' @export
#'
#' @examples
#' counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(normcounts = counts))
#' sce$pseudotime = seq_len(10)
#' as.AtgnatData(sce, pseudotime_slot="pseudotime", n_bins=3)
# TODO make bins and genes hidden with . and then add setters/getters?
AtgnatData = setClass(
  Class = "AtgnatData",
  slots = list(
    pseudobulks = "data.frame",
    bins = "numeric",
    genes = "character"
  )
)

#' @title Show an AtgnatData object
#'
#' @concept atgnat-object
#'
#' @param object an [AtgnatData] object
#' @export
#' @inherit AtgnatData-class examples
setMethod(f = "show",
          signature = "AtgnatData",
          definition = function(object){

            output = c(
              "Atgnat Data with:",
              paste("\tpseudobulks count:", ncol(object@pseudobulks)),
              paste("\ttotal genes:", nrow(object@pseudobulks)),
              paste("\tbins:", list(object@bins)),
              paste("\tselected genes:", list(object@genes), "\n")
            )

            cat(paste(output, collapse = '\n'))
          })

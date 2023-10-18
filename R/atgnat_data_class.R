#' Atgnat Data Object
#'
#' For creation details, see [as.AtgnatData()]
#'
#' @concept atgnat-object
#'
#' @slot pseudobulk_bins list of [data.frame]s. Each item is a normalised count
#' matrix representing a bin, where a column is a cell in the bin and each row is a gene.
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
    pseudobulk_bins = "list",
    bins = "numeric",
    genes = "character"
  )
)

# TODO validation of pseudobulk bins being a list of matrices

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
              paste("\tbins:", list(object@bins)),
              paste("\tselected genes:", list(object@genes), "\n")
            )

            cat(paste(output, collapse = '\n'))
          })

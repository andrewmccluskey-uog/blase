#' Blase Data Object
#'
#' For creation details, see [as.BlaseData()]
#'
#' @concept blase-object
#'
#' @slot pseudobulk_bins list of [data.frame]s. Each item is a normalised
#' count matrix representing a bin, where a column is a cell in the bin
#' and each row is a gene.
#' @slot bins list. A list of bin names for each timepoint.
#' @slot genes list. A list of the genes selected for
#' discriminating timepoints.
#'
#' @return An [BlaseData] object
#' @export
#'
#' @examples
#' counts <- matrix(rpois(100, lambda = 10), ncol = 10, nrow = 10)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(normcounts = counts)
#' )
#' sce$pseudotime <- seq_len(10)
#' data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 3)
#' genes(data) <- as.character(seq_len(10))
#' 
#' genes(data)
BlaseData <- setClass(
    Class = "BlaseData",
    slots = list(
        pseudobulk_bins = "list",
        bins = "numeric",
        genes = "character"
    )
)

# TODO validation of pseudobulk bins being a list of matrices

#' @title Show an BlaseData object
#'
#' @concept blase-object
#'
#' @param object an [BlaseData] object
#' @export
#' @inherit BlaseData-class examples
setMethod(
    f = "show",
    signature = "BlaseData",
    definition = function(object) {
        output <- c(
            "Blase Data with:",
            paste("\tbins:", list(object@bins)),
            paste("\tselected genes:", list(object@genes), "\n")
        )

        cat(paste(output, collapse = "\n"))
    }
)

#' @title Get genes of a BLASE Data object.
#'
#' @concept blase-object
#'
#' @rdname genes-getter
#' @param x a [BlaseData] object
#' @export
#' @inherit BlaseData-class examples
setGeneric("genes", function(x) standardGeneric("genes"))

#' @rdname genes-getter
setMethod(
  f = "genes", 
  signature = "BlaseData", 
  definition = function(x) x@genes)

#' @title Set genes of a BLASE Data object.
#'
#' @concept blase-object
#'
#' @rdname genes-setter
#' @param x a [BlaseData] object
#' @param value new value for genes slot, should be a vector of strings
#' @export
#' @inherit BlaseData-class examples
setGeneric("genes<-", function(x, value) standardGeneric("genes<-"))

#' @rdname genes-setter
setMethod(
  f = "genes<-", 
  signature = "BlaseData", 
  definition = function(x, value) {
    x@genes <- value
    validObject(x)
    x
  })


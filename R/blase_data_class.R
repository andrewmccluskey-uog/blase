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
#' @return A [BlaseData] object
#' @export
#'
#' @examples
#' counts <- matrix(rpois(100, lambda = 10), ncol = 10, nrow = 10)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(normcounts = counts)
#' )
#' sce$pseudotime <- seq_len(10) - 1
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


#' @title Show an BlaseData object
#'
#' @concept blase-object
#' @rdname show_blase_object
#'
#' @param object a [BlaseData] object
#'
#' @returns A character vector describing the BLASE object
#'
#' @importFrom methods show
#'
#' @export
#' @inherit BlaseData-class examples
setMethod(
    f = "show",
    signature = "BlaseData",
    definition = function(object) {
        output <- c(
            "Blase Data with:",
            paste0("\tbins: ", list(bins(object))),
            paste0("\tselected genes: ", list(genes(object)), "\n")
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
#'
#' @returns The vector of genes a BLASE object will use for mappings.
#'
#' @export
#' @inherit BlaseData-class examples
setGeneric("genes", function(x) standardGeneric("genes"))

#' @rdname genes-getter
setMethod(
    f = "genes",
    signature = "BlaseData",
    definition = function(x) x@genes
)

#' @title Set genes of a BLASE Data object.
#'
#' @concept blase-object
#'
#' @rdname genes-setter
#' @param x a [BlaseData] object
#' @param value Vector of strings. The new value for genes slot
#'
#' @returns Nothing
#'
#' @importFrom methods validObject
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
    }
)

#' @title Get bins of a BLASE Data object.
#' @keywords internal
#'
#' @concept blase-object
#'
#' @rdname bins-getter
#' @param x a [BlaseData] object
#'
#' @returns vector of integers. The bin names in the
#' BLASE Data object.
#'
#' @inherit BlaseData-class examples
setGeneric("bins", function(x) standardGeneric("bins"))

#' @rdname bins-getter
setMethod(
  f = "bins",
  signature = "BlaseData",
  definition = function(x) x@bins
)

#' @title Set genes of a BLASE Data object.
#' @keywords internal
#'
#' @concept blase-object
#'
#' @rdname bins-setter
#' @param x a [BlaseData] object
#' @param value Vector of integers The new value for bins slot
#'
#' @returns Nothing
#'
#' @importFrom methods validObject
setGeneric("bins<-", function(x, value) standardGeneric("bins<-"))

#' @rdname bins-setter
setMethod(
  f = "bins<-",
  signature = "BlaseData",
  definition = function(x, value) {
    x@bins <- value
    validObject(x)
    x
  }
)

#' @title Get pseudobulk bins of a BLASE Data object.
#' @keywords internal
#'
#' @concept blase-object
#'
#' @rdname pseudobulk-bins-getter
#' @param x a [BlaseData] object
#'
#' @returns List of dataframes. Each dataframe is the normalised counts of
#' cells by genes in each pseudotime bin. List index is the pseudotime bin.
#'
#' @inherit BlaseData-class examples
setGeneric("pseudobulk_bins", function(x) standardGeneric("pseudobulk_bins"))

#' @rdname pseudobulk-bins-getter
setMethod(
  f = "pseudobulk_bins",
  signature = "BlaseData",
  definition = function(x) x@pseudobulk_bins
)

#' @title Set genes of a BLASE Data object.
#' @keywords internal
#'
#' @concept blase-object
#'
#' @rdname pseudobulk-bins-setter
#' @param x a [BlaseData] object
#' @param value List of dataframes. Each dataframe is the normalised counts of
#' cells by genes in each pseudotime bin. List index is the pseudotime bin.
#'
#' @returns Nothing
#'
#' @importFrom methods validObject
setGeneric("pseudobulk_bins<-", function(x, value) {
  standardGeneric("pseudobulk_bins<-")
})

#' @rdname pseudobulk-bins-setter
setMethod(
  f = "pseudobulk_bins<-",
  signature = "BlaseData",
  definition = function(x, value) {
    x@pseudobulk_bins <- value
    validObject(x)
    x
  }
)


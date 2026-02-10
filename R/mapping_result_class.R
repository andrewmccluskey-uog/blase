#' Blase Mapping Result
#'
#' Created by [map_best_bin()]
#'
#' @concept mapping-result-object
#'
#' @param bulk_name String. The name of the bulk sample being mapped.
#' @param best_bin Integer. The bin that best matched the bulk sample.
#' @param best_correlation Decimal. The spearman's rho that the test geneset had
#' between the winning bin and the bulk.
#' @param top_2_distance Decimal. The absolute difference between the best and
#' second best mapping buckets. Higher indicates a less doubtful mapping.
#' @param strong_mapping Boolean. TRUE when the mapped bin's lower
#' bound is higher than the maximum upper bound of the other bins.
#' @param history A dataframe of the correlation score (decimal) and
#' confidence bounds (decimal pairs) for each bin.
#' Access with `mapping_history()`
#' @param bootstrap_iterations Integer. The number of iterations used during
#' the bootstrap.
#'
#' @return A MappingResult object
#'
#' @importFrom methods new
#'
#' @seealso [map_best_bin()]
#'
#' @examples
#' counts_matrix <- matrix(
#'     c(seq_len(120) / 10, seq_len(120) / 5),
#'     ncol = 48, nrow = 5
#' )
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(
#'     normcounts = counts_matrix, logcounts = log(counts_matrix)
#' ))
#' colnames(sce) <- seq_len(48)
#' rownames(sce) <- as.character(seq_len(5))
#' sce$cell_type <- c(rep("celltype_1", 24), rep("celltype_2", 24))
#'
#' sce$pseudotime <- seq_len(48) - 1
#' blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 4)
#' genes(blase_data) <- as.character(seq_len(5))
#'
#' bulk_counts <- matrix(seq_len(15) * 10, ncol = 3, nrow = 5)
#' colnames(bulk_counts) <- c("A", "B", "C")
#' rownames(bulk_counts) <- as.character(seq_len(5))
#'
#' # Map to bin
#' result <- map_best_bin(blase_data, "B", bulk_counts)
#' result
#'
#' # Map all bulks to bin
#' results <- map_all_best_bins(blase_data, bulk_counts)
#'
#' # Plot Heatmap
#' plot_mapping_result_heatmap(list(result))
#'
#' # Plot Correlation
#' plot_mapping_result_corr(result)
#'
#' # Plot populations
#' sce <- assign_pseudotime_bins(
#'     sce,
#'     pseudotime_slot = "pseudotime", n_bins = 4
#' )
#' plot_bin_population(sce, best_bin(result), group_by_slot = "cell_type")
#'
#' # Getters
#' bulk_name(result)
#' best_bin(result)
#' best_correlation(result)
#' top_2_distance(result)
#' strong_mapping(result)
#' mapping_history(result)
#' bootstrap_iterations(result)
#'
#' # Setters
#' bulk_name(result) <- "New Name"
MappingResult <- function(
    bulk_name,
    best_bin,
    best_correlation,
    top_2_distance,
    strong_mapping,
    history,
    bootstrap_iterations,
    metric="spearman") {

  return(methods::new("MappingResult",
          bulk_name = bulk_name,
          best_bin = best_bin,
          best_correlation = best_correlation,
          top_2_distance = top_2_distance,
          strong_mapping = strong_mapping,
          history = history,
          bootstrap_iterations = bootstrap_iterations,
          metric = metric
  ))

}

#' @export
#' @rdname MappingResult
setClass(
    Class = "MappingResult",
    slots = list(
        bulk_name = "ANY",
        best_bin = "numeric",
        best_correlation = "numeric",
        top_2_distance = "numeric",
        strong_mapping = "logical",
        history = "data.frame",
        bootstrap_iterations = "numeric",
        metric = "character"
    )
)

#' @title Show an MappingResult object
#'
#' @concept mapping-result-object
#'
#' @param object an [MappingResult] object
#'
#' @returns A character vector describing the Mapping Result object
#'
#' @importFrom methods show
#' @export
#' @inherit MappingResult examples
setMethod(
    f = "show",
    signature = "MappingResult",
    definition = function(object) {
        non_top_mapping_best_upper_bound <- max(
            mapping_history(object)[
              mapping_history(object)$bin != best_bin(object), ]$upper_bound)

        output <- c(
            paste0(
                "MappingResult for '",
                bulk_name(object),
                "': best_bin=", best_bin(object),
                " correlation=", best_correlation(object),
                " top_2_distance=", top_2_distance(object)
            ),
            paste(
                "\t Strong Result:",
                strong_mapping(object),
                "(next max upper ",
                non_top_mapping_best_upper_bound,
                ")"
            ),
            paste(
                "\t with history for scores against",
                nrow(mapping_history(object)),
                " bins"
            ),
            paste(
                "\t Bootstrapped with",
                bootstrap_iterations(object),
                "iterations\n"
            )
        )

        cat(paste(output, collapse = "\n"))
    }
)

#' @title Get name of bulk of a BLASE Mapping Results object.
#'
#' @concept mapping-result-object
#'
#' @rdname mapping-result-bulk-name-getter
#' @param x a [MappingResult] object
#' @returns String. The name of the bulk used to map against.
#' @export
#' @inherit MappingResult examples
setGeneric("bulk_name", function(x) standardGeneric("bulk_name"))

#' @rdname mapping-result-bulk-name-getter
setMethod(
    f = "bulk_name",
    signature = "MappingResult",
    definition = function(x) x@bulk_name
)

#' @title Set name of bulk of a BLASE Mapping Results object.
#'
#' @concept mapping-result-object
#'
#' @rdname bulk_name-setter
#' @param x a [MappingResult] object
#' @param value String. The name of the bulk used to map against.
#'
#' @returns Nothing
#'
#' @importFrom methods validObject
#' @export
#' @inherit BlaseData-class examples
setGeneric("bulk_name<-", function(x, value) standardGeneric("bulk_name<-"))

#' @rdname bulk_name-setter
setMethod(
  f = "bulk_name<-",
  signature = "MappingResult",
  definition = function(x, value) {
    x@bulk_name <- value
    validObject(x)
    x
  }
)


#' @title Get best bin of a BLASE Mapping Results object.
#'
#' @concept mapping-result-object
#'
#' @rdname mapping-result-best-bin-getter
#' @param x a [MappingResult] object
#' @returns Integer. The best bin ID of this mapping
#' @export
#' @inherit MappingResult examples
setGeneric("best_bin", function(x) standardGeneric("best_bin"))

#' @rdname mapping-result-best-bin-getter
setMethod(
    f = "best_bin",
    signature = "MappingResult",
    definition = function(x) x@best_bin
)

#' @title Get best correlation of a BLASE Mapping Results object.
#'
#' @concept mapping-result-object
#'
#' @rdname mapping-result-best-correlation-getter
#' @param x a [MappingResult] object
#' @returns Decimal. The highest correlation value of this mapping
#' @export
#' @inherit MappingResult examples
setGeneric("best_correlation", function(x) standardGeneric("best_correlation"))

#' @rdname mapping-result-best-correlation-getter
setMethod(
    f = "best_correlation",
    signature = "MappingResult",
    definition = function(x) x@best_correlation
)

#' @title Get the difference in correlation between the top 2 most correlated
#' bins for a BLASE Mapping Results object.
#'
#' @concept mapping-result-object
#'
#' @rdname mapping-result-top-2-distance-getter
#' @param x a [MappingResult] object
#' @returns Decimal. The difference in correlation between
#' the top 2 most correlated bins for this mapping.
#' @export
#' @inherit MappingResult examples
setGeneric(
    "top_2_distance",
    function(x) standardGeneric("top_2_distance")
)

#' @rdname mapping-result-top-2-distance-getter
setMethod(
    f = "top_2_distance",
    signature = "MappingResult",
    definition = function(x) x@top_2_distance
)

#' @title Get if the result is strong for a BLASE Mapping Results object.
#'
#' @concept mapping-result-object
#'
#' @rdname mapping-result-strong-mapping-getter
#' @param x a [MappingResult] object
#' @returns Boolean. TRUE if the result is strong, otherwise FALSE
#' @export
#' @inherit MappingResult examples
setGeneric(
    "strong_mapping",
    function(x) standardGeneric("strong_mapping")
)

#' @rdname mapping-result-strong-mapping-getter
setMethod(
    f = "strong_mapping",
    signature = "MappingResult",
    definition = function(x) x@strong_mapping
)

#' @title Get the mapping history for a BLASE Mapping Results object.
#'
#' @concept mapping-result-object
#'
#' @rdname mapping-result-history-getter
#' @param x a [MappingResult] object
#' @returns The mapping history of this mapping, in a data frame.
#' @export
#' @inherit MappingResult examples
setGeneric("mapping_history", function(x) standardGeneric("mapping_history"))

#' @rdname mapping-result-history-getter
setMethod(
    f = "mapping_history",
    signature = "MappingResult",
    definition = function(x) x@history
)

#' @title Get the number of bootstrap iterations
#' performed for a BLASE Mapping Results object.
#'
#' @concept mapping-result-object
#'
#' @rdname mapping-result-bootstrap-iterations-getter
#' @param x a [MappingResult] object
#' @returns Integer. The number of iterations performed for this mapping.
#' @export
#' @inherit MappingResult examples
setGeneric(
    "bootstrap_iterations",
    function(x) standardGeneric("bootstrap_iterations")
)

#' @rdname mapping-result-bootstrap-iterations-getter
setMethod(
    f = "bootstrap_iterations",
    signature = "MappingResult",
    definition = function(x) x@bootstrap_iterations
)

#' @title Get the mapping history for a BLASE Mapping Results object.
#'
#' @concept mapping-result-object
#'
#' @rdname mapping-result-metric-getter
#' @param x a [MappingResult] object
#' @returns a String, the metric used to calculate the result.
#' @export
#' @inherit MappingResult examples
setGeneric("metric", function(x) standardGeneric("metric"))

#' @rdname mapping-result-metric-getter
setMethod(
  f = "metric",
  signature = "MappingResult",
  definition = function(x) x@metric
)

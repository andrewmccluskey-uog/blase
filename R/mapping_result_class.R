#' Atgnat Mapping Result
#'
#' Created by [map_best_bin()]
#'
#' @concept mapping
#'
#' @slot best_bin The bin that best matched the bulk sample.
#' @slot best_correlation The spearman's rho that the test geneset had between the winning bin and the bulk.
#' @slot top_2_distance The absolute difference between the best and second best mapping buckets. Higher indicates a less doubtful mapping.
#' @slot history A dataframe of the correlation score for each bin.
#'
#' @return A [MappingResult] object
#' @export
#'
#' @seealso [map_best_bin()]
#'
#' @examples
#'
MappingResult = setClass(
  Class = "MappingResult",
  slots = list(
    bulk_name = "ANY",
    best_bin = "numeric",
    best_correlation = "numeric",
    top_2_distance = "numeric",
    history = "data.frame"
  )
)

#' @title Show an MappingResult object
#'
#' @concept mapping
#'
#' @param object an [MappingResult] object
#' @export
setMethod(f = "show",
          signature = "MappingResult",
          definition = function(object){

            output = c(
              paste0("MappingResult for '", object@bulk_name, "':",
                    " best_bin=", object@best_bin,
                     " correlation=", object@best_correlation,
                     " top_2_distance=", object@top_2_distance),
              paste("\t with history for scores against", nrow(object@history), " bins\n")
            )

            cat(paste(output, collapse = '\n'))
          })

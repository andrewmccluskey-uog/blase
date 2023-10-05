#' @title Plot a mapping result
#' @rdname plot_mapping_result
#' @param x An object to plot on.
#' @param y The [MappingResult] object to plot
#' @param ... additional arguments passed to object-specific methods.
#'
#' @export
setGeneric(name = "plot_mapping_result",
           signature = c(x="x", y="y"),
           def = function(x, y, ...) standardGeneric("plot_mapping_result"))

#' @rdname plot_mapping_result
#'
#' @export
setMethod(
  f = "plot_mapping_result",
  signature = c(x="SingleCellExperiment", y="MappingResult"),
  definition = function(x, y){
    # TODO the plotting
  }
)

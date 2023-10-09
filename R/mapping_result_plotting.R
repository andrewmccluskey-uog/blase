#' @title Plot a mapping result
#'
#' @concept mapping
#'
#' @rdname plot_mapping_result
#' @param x An object to plot on.
#' @param y The [MappingResult] object to plot
#' @param ... additional arguments passed to object-specific methods.
#'
#'
#' @export
setGeneric(name = "plot_mapping_result",
           signature = c(x="x", y="y"),
           def = function(x, y, ...) standardGeneric("plot_mapping_result"))

#' @rdname plot_mapping_result
#'
#' @param group_by_slot The slot in the [SingleCellExperiment::SingleCellExperiment]
#' to be used as the coloring for the output plot. Passed to [scater::plotUMAP()] as
#' `colour_by`, and will be used to produce a bar chart of populations in the
#' best mapped bin.
#'
#' @import scater
#'
#' @export
setMethod(
  f = "plot_mapping_result",
  signature = c(x="SingleCellExperiment", y="MappingResult"),
  definition = function(x, y, group_by_slot){

    bin_sym = ggplot2::sym("bin")
    correlation_sym = ggplot2::sym("correlation")
    var1_sym = ggplot2::sym("Var1")
    freq_sym = ggplot2::sym("Freq")

    best_bin_population_data = as.data.frame( table( x[,x$pseudotime_bin==y@best_bin]@colData[[group_by_slot]] ))

    # TODO offer this with more than just umaps
    gridExtra::grid.arrange(
      scater::plotUMAP(x, colour_by="pseudotime_bin"),
      scater::plotUMAP(x, colour_by=group_by_slot),
      scater::plotUMAP(x[,x$pseudotime_bin==y@best_bin], colour_by="pseudotime_bin"),
      scater::plotUMAP(x[,x$pseudotime_bin==y@best_bin], colour_by=group_by_slot),
      ggplot2::ggplot(y@history, ggplot2::aes(x={{bin_sym}}, y={{correlation_sym}})) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept=y@best_correlation, linetype="dashed") +
        ggplot2::geom_vline(xintercept=y@best_bin, linetype="dashed"),
      ggplot2::ggplot(best_bin_population_data[best_bin_population_data$Freq>0,], ggplot2::aes(x={{var1_sym}}, y={{freq_sym}})) + ggplot2::geom_bar(stat="identity"),
      ncol=2,
      top = grid::textGrob(paste0(y@bulk_name, ": Bin ", y@best_bin, ", Cor ", round(y@best_correlation, 4),", distance ", y@top_2_distance), gp=grid::gpar(fontsize=20,font=3))
    )
  }
)

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
  definition = function(x, y, group_by_slot){

    best_bin_population_data = as.data.frame( table( subset(x, , pseudotime_bin==y@best_bin)@colData[[group_by_slot]] ) )

    # TODO offer this with more than just umaps
    gridExtra::grid.arrange(
      scater::plotUMAP(x, text_by="pseudotime_bin", colour_by="pseudotime_bin"),
      scater::plotUMAP(x, colour_by=group_by_slot),
      scater::plotUMAP(x[,x$pseudotime_bin==y@best_bin], colour_by="pseudotime_bin"),
      scater::plotUMAP(x[,x$pseudotime_bin==y@best_bin], colour_by=group_by_slot),
      ggplot2::ggplot(y@history, aes(x=bin, y=correlation)) +
        geom_line() +
        geom_hline(yintercept=y@best_correlation, linetype="dashed") +
        geom_vline(xintercept=y@best_bin, linetype="dashed"),
      ggplot2::ggplot(best_bin_population_data[best_bin_population_data$Freq>0,], aes(x=Var1, y=Freq)) + geom_bar(stat="identity"),
      ncol=2,
      top = grid::textGrob(paste0(y@bulk_name, ": Bin ", y@best_bin, ", Cor ", round(y@best_correlation, 4),", distance ", y@top_2_distance), gp=grid::gpar(fontsize=20,font=3))
    )
  }
)

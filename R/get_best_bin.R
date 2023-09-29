#' Map the best matching SC bin for a bulk sample
#'
#' @param pseudotime_sce The SingleCellObject holding the bins.
#' @param pseudotime_gene_ratios The dataframe of gene ratios per bin.
#' @param bulk_id The sample id of the bulk to analyse
#' @param bulk_data The whole bulk read matrix
#' @param make_plot Whether or not to draw a plot of the results.
#'
#' @return A dataframe with one row, containing results from the mapping process.
#' * Bin: the bin that best matched the bulk sample.
#' * Corr: the spearman's rho that the test geneset had between the winning bin and the bulk.
#' * top_2_distance: the absolute difference between the best and second best mapping buckets. Higher indicates a less doubtful mapping.
#' * history: a dataframe of the correlation score for each bin.
#' @export
#'
#' @examples
map_best_bin <- function(pseudotime_sce, pseudotime_gene_ratios, bulk_id, bulk_data, make_plot=FALSE) {

  sum_for_top_genes = bulk_data[colnames(pseudotime_gene_ratios),bulk_id]
  best_cor = -1
  best_i = 0
  correlations_history = data.frame()
  for (i in seq_len(length(rownames(pseudotime_gene_ratios))) ) {
    bin_ratios = t(pseudotime_gene_ratios)[,i]
    corr <- stats::cor.test(bin_ratios, sum_for_top_genes, method = 'spearman')
    if (corr$estimate > best_cor) {
      best_cor = unname(corr$estimate)
      best_i = i
    }
    correlations_history <- rbind(correlations_history, c(i, unname(corr$estimate)))
    i=i+1
  }
  colnames(correlations_history) = c('bin', 'correlation')
  best_bin_population_data = as.data.frame(table(subset(pseudotime_sce, , pseudotime_sce@colData[["pseudotime_bin"]]==best_i)$Stages))

  top2 = utils::head(sort(correlations_history$correlation, decreasing=TRUE),n=2)
  # TODO round this to 4 decimal places not significant figures
  distance_between_top_2_corrs = signif(top2[1]-top2[2],2)

  # TODO resolve build issues with scater
  # TODO split out into different function
  #if (make_plot == TRUE) {
  #  gridExtra::grid.arrange(
  #    scater::plotUMAP(pseudotime_sce, text_by="pseudotime_bin", colour_by="pseudotime_bin"),
  #    scater::plotUMAP(pseudotime_sce, colour_by="Stages"),
  #    scater::plotUMAP(pseudotime_sce[,pseudotime_sce$pseudotime_bin==best_i], colour_by="pseudotime_bin"),
  #    scater::plotUMAP(pseudotime_sce[,pseudotime_sce$pseudotime_bin==best_i], colour_by="Stages"),
  #    ggplot2::ggplot(correlations_history, aes(x={{ggplot2::sym("bin")}}, y={{ggplot2::sym("correlation")}})) +
  #      ggplot2::geom_line() +
  #      ggplot2::geom_hline(yintercept=best_cor, linetype="dashed") +
  #      ggplot2::geom_vline(xintercept=best_i, linetype="dashed"),
  #    ggplot2::ggplot(best_bin_population_data[best_bin_population_data$Freq>0,], aes(x={{ggplot2::sym("Var1")}}, y={{ggplot2::sym("Freq")}})) + geom_bar(stat="identity"),
  #    ncol=2,
  #    top = gridExtra::textGrob(paste(bulk_sample, "( Bin", best_i, ", Cor", signif(best_cor, 2),", distance", distance_between_top_2_corrs, ")"), gp=gpar(fontsize=20,font=3))
  #  )
  #}

  result = data.frame(bin=best_i, correlation=best_cor, top_2_distance=distance_between_top_2_corrs, history=correlations_history)
  return(result)
}

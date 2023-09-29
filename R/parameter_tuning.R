#' Identify the Best Parameters For Your Dataset
#'
#' @param sce The Single Cell Experiment object with the counts
#' @param associationTestResult The results from your tradeseq associationTest
#' @param lineage The lineage to use
#' @param bins_count_range The n_bins list to try out
#' @param gene_count_range The n_genes list to try out
#' @param split_by How to do the pseudotime bucket creation
#'
#' @return A dataframe of the results.
#' * bin_count: The bin count for this attempt
#' * gene_count: The top n genes to use for this attempt
#' * worst_specificity: The worst specificity for these parameters
#' * mean_specificity: The mean specificity for these parameters
#'
#' @seealso [plot_find_best_params_results()]
#'
#' @export
#'
#' @examples
# TODO: consider allowing this to take just a list of genes not an association test
#   result to make it a bit more reusable. If we ask it to be ordered, we can
#   just take the top n.
find_best_params <- function(sce, associationTestResult, lineage=NA, bins_count_range=c(5,10,15,20,25,30,35), gene_count_range=c(20,30,40,45,50,55,60,70,80), split_by="pseudotime_range") {

  slingshot_pseudotime_slot = "slingPseudotime_1"

  if (!is.na(lineage)) {
    slingshot_pseudotime_slot = paste0("slingPseudotime_", lineage)
  }

  results = data.frame(gene_count=c(), bin_count=c(), worst_specificity=c(), mean_specificity=c())
  # TODO consider parallelising this
  # https://www.bioconductor.org/packages/devel/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.html#parallel-looping-vectorized-and-aggregate-operations
  # https://cran.r-project.org/web/packages/foreach/index.html
  for (bin_count in bins_count_range) {
    print(bin_count)
    sce = create_pseudotime_bins(sce, bin_count, pseudotime_slot = slingshot_pseudotime_slot,  split_by = split_by)

    for (genes_count in gene_count_range) {
      gene_list = get_top_n_genes(associationTestResult, n_genes=genes_count, lineage=lineage)
      pseudotime_bins_ratios_of_top_n_genes_df = get_pseudotime_gene_ratios(sce, gene_list)
      res = evaluate_parameters(sce, pseudotime_bins_ratios_of_top_n_genes_df, make_plot=FALSE)
      results = rbind(results, data.frame(bin_count=c(bin_count), gene_count=c(genes_count), worst_specificity=c(res[1]), mean_specificity=c(res[2])))
    }
  }

  return(results)

}

#' Plot the results of the search for good parameters
#'
#' @param find_best_params_results Results dataframe from [find_best_params()]
#' @param bin_count_colors Optional, custom bin count color scheme.
#' @param gene_count_colors Optional, custom gene count color scheme.
#'
#' @seealso [find_best_params()]
#'
#' @export
#'
#' @examples
#' best_param_results = data.frame(
#'   gene_count=c(10,10,20,20),
#'   bin_count=c(15,25,15,25),
#'   worst_specificity=c(0.01,0.02,0.015,0.025),
#'   mean_specificity=c(0.2,0.1,0.11,0.23))
#' plot_find_best_params_results(best_param_results)
plot_find_best_params_results <- function(find_best_params_results, bin_count_colors=viridis::scale_color_viridis(option="viridis"), gene_count_colors=viridis::scale_color_viridis(option="magma")) {
  gene_count = ggplot2::sym("gene_count")
  bin_count = ggplot2::sym("bin_count")
  worst_specificity = ggplot2::sym("worst_specificity")
  mean_specificity = ggplot2::sym("mean_specificity")

  gridExtra::grid.arrange(
    ggplot2::ggplot(find_best_params_results, ggplot2::aes(x={{gene_count}}, y={{worst_specificity}}, color={{bin_count}})) + ggplot2::geom_point() + bin_count_colors,
    ggplot2::ggplot(find_best_params_results, ggplot2::aes(x={{bin_count}}, y={{worst_specificity}}, color={{gene_count}})) + ggplot2::geom_point() + gene_count_colors,
    ggplot2::ggplot(find_best_params_results, ggplot2::aes(x={{gene_count}}, y={{mean_specificity}}, color={{bin_count}})) + ggplot2::geom_point() + bin_count_colors,
    ggplot2::ggplot(find_best_params_results, ggplot2::aes(x={{bin_count}}, y={{mean_specificity}}, color={{gene_count}})) + ggplot2::geom_point() + gene_count_colors,
    ncol=2)
}

#' Evaluate n_bins and n_genes for bin mapping
#'
#' @description Will use the n_bins and n_genes implied by the `sce` and
#' `pseudotime_bins_top_n_genes_df` parameters and return quality metrics and
#' an optional chart.
#'
#' @param sce The Single Cell Experiment object to use.
#' @param pseudotime_bins_top_n_genes_df The dataframe with the expression for
#' the top genes per bin.
#' @param bins_slot The slot that the bins are in.
#' @param make_plot Whether or not to render the plot showing the correlations
#' for each pseudobulk bin when we try to map the given bin.
#' @param plot_columns How many columns to use in the plot.
#'
#' @return A vector of length 2:
#' * "worst top 2 distance" containing the lowest difference between the absolute values of the
#'  top 2 most correlated bins for each bin. Higher is better for differentiating.
#' * "mean top 2 distance" containing the mean top 2 distance across the entire set of genes and bins.
#'  Higher is better for differentiation, but it should matter less than the worst value.
#' @export
#'
#' @examples
evaluate_parameters <- function(sce, pseudotime_bins_top_n_genes_df, bins_slot="pseudotime_bin", make_plot=FALSE, plot_columns=4) {

  bin_ids = as.numeric(names(table(sce@colData[[bins_slot]])))

  pseudobulks = PRIVATE_get_raw_pseudobulks(sce, bin_ids)

  results.best_bin = c()
  results.best_corr = c()
  results.history = c()
  results.specificity = c()

  for (i in bin_ids) {
    res = map_best_bin(sce, pseudotime_bins_top_n_genes_df, i, pseudobulks, make_plot=FALSE)
    results.best_bin = append(results.best_bin, c(res[[1,1]]))
    results.best_corr = append(results.best_corr, c(res[[1,2]]))
    results.specificity = append(results.specificity, c(res[[1,3]]))
    results.history = append(results.history, c(res[,3:4]))
  }

  worst_specificity = min(results.specificity)
  mean_specificity = mean(results.specificity)

  if (make_plot == TRUE) {
    plots=list()

    for (i in bin_ids) {
      plots[[i]] = PRIVATE_plot_history(i, results.best_bin, results.best_corr, results.history, results.specificity)
    }

    gridExtra::grid.arrange(
      top=grid::textGrob(paste(length(colnames(pseudotime_bins_top_n_genes_df)), "genes and worst specificity:", signif(worst_specificity, 2)),gp=grid::gpar(fontsize=20,font=3)),
      grobs=plots,
      ncol=plot_columns
    )
  }

  return(c(worst_specificity, mean_specificity))

}

#' Evaluate Top Genes
#'
#' Shows plots over bins of expression of the top n genes. This is designed to help
#' identify if you have selected genes that vary over the pseudotime you have chosen
#' bins to exist over. Uses the logcounts of the SCE.
#'
#' @param sce The Single Cell Experiment to get bins and expression from.
#' @param pseudotime_bins_top_n_genes_df The dataframe of genes expression by bin.
#' @param bins_slot The slot to look for the bin IDs in.
#' @param n_genes_to_plot The number of genes to plot.
#' @param plot_columns The number of columns to plot the grid with. Best as a
#' divisor of `n_genes_to_plot`.
#'
#' @export
#'
#' @examples
evaluate_top_n_genes <- function(sce, pseudotime_bins_top_n_genes_df, bins_slot="pseudotime_bin", n_genes_to_plot=16, plot_columns=4) {

  bin_ids = as.numeric(names(table(sce@colData[bins_slot])))
  pseudobulks = PRIVATE_get_raw_pseudobulks(sce, bin_ids)

  plots = list()
  for (i in seq_len(n_genes_to_plot)) {
    plots[[i]] = PRIVATE_plot_gene_over_bins(pseudobulks, colnames(pseudotime_bins_top_n_genes_df)[i])
  }

  gridExtra::grid.arrange(
    top=grid::textGrob(paste(ncol(pseudotime_bins_top_n_genes_df), "genes"),gp=grid::gpar(fontsize=20,font=3)),
    grobs=plots,
    ncol=plot_columns
  )

}

PRIVATE_get_raw_pseudobulks = function(sce, bin_ids) {
  pseudobulks = as.data.frame(Matrix::rowSums(SingleCellExperiment::logcounts(sce)))
  for (i in bin_ids) {

    pb = as.data.frame(Matrix::rowSums(SingleCellExperiment::logcounts(subset(sce, , sce@colData[["pseudotime_bin"]] == i))))
    colnames(pb) = i
    pseudobulks = cbind(pseudobulks, pb)
  }
  pseudobulks = pseudobulks[2:length(pseudobulks)]
  return(pseudobulks)
}

PRIVATE_plot_history = function(i, bin, corr, history, specificity){
  return(ggplot2::ggplot(as.data.frame(history[(i*2-1):(i*2)]), ggplot2::aes(x={{ggplot2::sym("bin")}}, y={{ggplot2::sym("correlation")}})) +
           ggplot2::ylim(-1,1) +
           ggplot2::ggtitle(paste0(bin[i], " (", signif(corr[i], 2) ,",",signif(specificity[i], 2),")")) +
           ggplot2::geom_line() +
           ggplot2::geom_hline(yintercept=corr[i], linetype="dashed") +
           ggplot2::geom_vline(xintercept=bin[i], linetype="dashed"))

}

PRIVATE_plot_gene_over_bins = function(pseudobulks, gene){

  expression =as.data.frame(t(pseudobulks[gene,]))
  colnames(expression) = "expr"
  expression$bin = seq_len(nrow(expression))

  return(ggplot2::ggplot(expression, ggplot2::aes(x={{ggplot2::sym("bin")}}, y={{ggplot2::sym("expr")}})) +
           ggplot2::ggtitle(gene) +
           ggplot2::geom_line())

}

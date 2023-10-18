#' Evaluate n_bins and n_genes for bin mapping
#'
#' @description Will use the n_bins and n_genes implied by the `sce` and
#' `pseudotime_bins_top_n_genes_df` parameters and return quality metrics and
#' an optional chart.
#'
#' @concept tuning
#'
#' @param atgnat_data The [AtgnatData] object to use.
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
#' @inherit find_best_params examples
evaluate_parameters <- function(atgnat_data, make_plot=FALSE, plot_columns=4) {

  bin_ids = atgnat_data@bins

  results.best_bin = c()
  results.best_corr = c()
  results.history = c()
  results.specificity = c()

  pseudobulked_bins = data.frame(lapply(seq_len(length(atgnat_data@pseudobulk_bins)), function (i) {
    x = atgnat_data@pseudobulk_bins[[i]]
    return(Matrix::rowMeans(x))
  }))
  colnames(pseudobulked_bins) = atgnat_data@bins

  for (i in bin_ids) {
    res = map_best_bin(atgnat_data, i, pseudobulked_bins, bootstrap_iterations = 0)
    results.best_bin = append(results.best_bin, c(res@best_bin))
    results.best_corr = append(results.best_corr, c(res@best_correlation))
    results.specificity = append(results.specificity, c(res@top_2_distance))
    results.history = append(results.history, c(res@history))
  }

  worst_specificity = min(results.specificity)
  mean_specificity = mean(results.specificity)

  if (make_plot == TRUE) {
    plots=list()

    for (i in bin_ids) {
      plots[[i]] = PRIVATE_plot_history(i, results.best_bin, results.best_corr, results.history, results.specificity)
    }

    gridExtra::grid.arrange(
      top=grid::textGrob(paste(length(atgnat_data@genes), "genes and worst specificity:", signif(worst_specificity, 2)),gp=grid::gpar(fontsize=20,font=3)),
      grobs=plots,
      ncol=plot_columns
    )
  }

  return(c(worst_specificity, mean_specificity))

}


#' Identify the Best Parameters For Your Dataset
#'
#' @concept tuning
#'
#' @param x The object to create `AtgnatData`` from
#' @param genelist The list of genes to use (ordered by descending goodness)
#' @param bins_count_range The n_bins list to try out
#' @param gene_count_range The n_genes list to try out
#' @param ... params to be passed to child functions, see [as.AtgnatData()]
#'
#' @return A dataframe of the results.
#' * bin_count: The bin count for this attempt
#' * gene_count: The top n genes to use for this attempt
#' * worst_specificity: The worst specificity for these parameters
#' * mean_specificity: The mean specificity for these parameters
#'
#' @seealso [plot_find_best_params_results()] for plotting the results of this function.
#'
#' @export
#'
#' @examples
#' ncells = 70
#' ngenes = 100
#' counts_matrix <- matrix(c(seq_len(3500)/10, seq_len(3500)/5), ncol=ncells, nrow=ngenes)
#' sce <- SingleCellExperiment::SingleCellExperiment(assays=list(
#'   normcounts=counts_matrix, logcounts=log(counts_matrix)))
#' colnames(sce) = seq_len(ncells)
#' rownames(sce) = as.character(seq_len(ngenes))
#' sce$cell_type = c(rep("celltype_1", ncells/2), rep("celltype_2", ncells/2))
#'
#' sce$pseudotime = seq_len(ncells)
#' genelist = as.character(seq_len(ngenes))
#'
#' # Finding the best params for the AtgnatData
#' best_params = find_best_params(sce, genelist, pseudotime_slot="pseudotime", split_by="cells")
#' best_params
#' plot_find_best_params_results(best_params)
#'
#' # Evaluating created AtgnatData
#' atgnat_data = as.AtgnatData(sce, pseudotime_slot="pseudotime", n_bins=10)
#' atgnat_data@genes = genelist[1:20]
#'
#' # Check specificity of parameters
#' evaluate_parameters(atgnat_data, make_plot = TRUE)
#' # Check gene expression over pseudotime
#' evaluate_top_n_genes(atgnat_data)
find_best_params <- function(x, genelist, bins_count_range=c(5,10,15,20,25,30,35), gene_count_range=c(20,30,40,45,50,55,60,70,80), ...) {

  results = data.frame(gene_count=c(), bin_count=c(), worst_specificity=c(), mean_specificity=c())
  # TODO consider parallelising this
  # https://www.bioconductor.org/packages/devel/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.html#parallel-looping-vectorized-and-aggregate-operations
  # https://cran.r-project.org/web/packages/foreach/index.html
  for (bin_count in bins_count_range) {
    atgnat_data = as.AtgnatData(x=x, n_bins=bin_count, ...)

    for (genes_count in gene_count_range) {
      atgnat_data@genes = genelist[1:genes_count]
      res = evaluate_parameters(atgnat_data, make_plot=FALSE)
      results = rbind(results, data.frame(bin_count=c(bin_count), gene_count=c(genes_count), worst_specificity=c(res[1]), mean_specificity=c(res[2])))
    }
  }

  return(results)

}

#' Plot the results of the search for good parameters
#'
#' @concept tuning
#'
#' @param find_best_params_results Results dataframe from [find_best_params()]
#' @param bin_count_colors Optional, custom bin count color scheme.
#' @param gene_count_colors Optional, custom gene count color scheme.
#'
#' @seealso [find_best_params()]
#'
#' @export
#'
#' @inherit find_best_params examples
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

#' Evaluate Top Genes
#'
#' Shows plots over bins of expression of the top n genes. This is designed to help
#' identify if you have selected genes that vary over the pseudotime you have chosen
#' bins to exist over. Uses the normcounts of the SCE.
#'
#' @concept tuning
#'
#' @param atgnat_data The [AtgnatData] to get bins and expression from.
#' @param n_genes_to_plot The number of genes to plot.
#' @param plot_columns The number of columns to plot the grid with. Best as a
#' divisor of `n_genes_to_plot`.
#'
#' @export
#'
#' @inherit find_best_params examples
evaluate_top_n_genes <- function(atgnat_data, n_genes_to_plot=16, plot_columns=4) {

  if (n_genes_to_plot > length(atgnat_data@genes)) {
    n_genes_to_plot = length(atgnat_data@genes)
  }

  pseudobulked_bins = data.frame(lapply(seq_len(length(atgnat_data@pseudobulk_bins)), function (i) {
    x = atgnat_data@pseudobulk_bins[[i]]
    return(Matrix::rowMeans(x))
  }))
  colnames(pseudobulked_bins) = atgnat_data@bins

  plots = list()
  for (i in seq_len(n_genes_to_plot)) {
    plots[[i]] = PRIVATE_plot_gene_over_bins(pseudobulked_bins, atgnat_data@genes[i])
  }

  gridExtra::grid.arrange(
    top=grid::textGrob(paste(length(atgnat_data@genes), "genes"),gp=grid::gpar(fontsize=20,font=3)),
    grobs=plots,
    ncol=plot_columns
  )

}

PRIVATE_plot_history = function(i, bin, corr, history, specificity){

  bin_sym = ggplot2::sym("bin")
  corr_sym = ggplot2::sym("correlation")

  return(ggplot2::ggplot(as.data.frame(history[(i*4-3):(i*4)]), ggplot2::aes(x={{bin_sym}}, y={{corr_sym}})) +
           ggplot2::ylim(-1,1) +
           ggplot2::ggtitle(paste0(bin[i], " (", signif(corr[i], 2) ,",",signif(specificity[i], 2),")")) +
           ggplot2::geom_line() +
           ggplot2::geom_hline(yintercept=corr[i], linetype="dashed") +
           ggplot2::geom_vline(xintercept=bin[i], linetype="dashed"))

}

PRIVATE_plot_gene_over_bins = function(pseudobulks, gene){

  bin_sym = ggplot2::sym("bin")
  expr_sym = ggplot2::sym("expr")

  expression =as.data.frame(t(pseudobulks[gene,]))
  colnames(expression) = "expr"
  expression$bin = seq_len(nrow(expression))

  return(ggplot2::ggplot(expression, ggplot2::aes(x={{bin_sym}}, y={{expr_sym}})) +
           ggplot2::ggtitle(gene) +
           ggplot2::geom_line())

}


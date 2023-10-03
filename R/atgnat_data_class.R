# TODO split methods into a) conversion.R b) mapping.R c) plotting.R d) tuning.R

#' Atgnat Data Object
#'
#' @slot pseudobulks data.frame. Each column is a timepoint sample and each row is a gene.
#' @slot bins list. A list of bin names for each timepoint.
#' @slot genes list. A list of the genes selected for discriminating timepoints.
#'
#' @return An [AtgnatData] object
#' @export
#'
#' @examples
#' library(SingleCellExperiment, quietly=TRUE)
#' library(atgnat)
#' counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(normcounts = counts))
#' sce$pseudotime = seq_len(10)
#' as.AtgnatData(sce, pseudotime_slot="pseudotime", n_bins=3)
# TODO make bins and genes hidden with . and then add setters/getters?
AtgnatData = setClass(
  Class = "AtgnatData",
  slots = list(
    pseudobulks = "data.frame",
    bins = "numeric",
    genes = "character"
  )
)

#' @title Show an AtgnatData object
#' @param object an [AtgnatData] object
#' @export
setMethod(f = "show",
          signature = "AtgnatData",
          definition = function(object){

            output = c(
              "Atgnat Data with:",
              paste("\tpseudobulks count:", ncol(object@pseudobulks)),
              paste("\ttotal genes:", nrow(object@pseudobulks)),
              paste("\tbins:", list(object@bins)),
              paste("\tselected genes:", object@genes)
            )

            cat(paste(output, collapse = '\n'))
          })

#' @title Conversion to AtgnatData
#' @rdname as.AtgnatData
#' @param x An object to take counts from
#' @param ... additional arguments passed to object-specific methods.
#'
#' @return An [AtgnatData] object
#'
#' @export
setGeneric(name = "as.AtgnatData",
           signature = c("x"),
           def = function(x, ...) standardGeneric("as.AtgnatData"))

#' @rdname as.AtgnatData
#'
#' @import methods
#'
#' @export
setMethod(
  f = "as.AtgnatData",
  signature = c(x="data.frame"),
  definition = function(x){
    return(methods::new("AtgnatData", pseudobulks = x, bins = colnames(x)))
  }
)


#' @rdname as.AtgnatData
#' @param pseudotime_slot The [SingleCellExperiment::SingleCellExperiment] slot containing the pseudotime values for each cell
#' @param n_bins The number of bins to create.
#' @param split_by The technique used to split the bins. The default `pseudotime_range` picks the
#' bin for a cell based on a constant range of pseudotime. `cells` picks the bin for a cell based on
#' an even number of cells per bin.
#'
#' @import methods
#'
#' @export
setMethod(
  f = "as.AtgnatData",
  signature = c(x="SingleCellExperiment"),
  definition = function(x, pseudotime_slot="slingPseudotime_1", n_bins=20, split_by="pseudotime_range"){

    if (is.na(match(split_by, c("pseudotime_range", "cells")))) {
      stop("split_by must be 'pseudotime_range' or 'cells'")
    }

    # TODO add check for pseudotime_slot existing

    pseudotime_sce = subset(x, , !is.na(x@colData[pseudotime_slot]))
    pseudotime = pseudotime_sce@colData[[pseudotime_slot]]

    c = SingleCellExperiment::normcounts(pseudotime_sce)

    if (split_by == "pseudotime_range") {
      min_pdt = 0
      max_pdt = ceiling(max(pseudotime))

      bin_size = max_pdt/n_bins
      bin_upper_limits = seq(bin_size, max_pdt, by=bin_size)

      ## Put cells into the right bins
      pseudotime_sce$pseudotime_bin = ceiling(pseudotime / bin_size)

      pseudotime_sce$pseudotime_bin[pseudotime_sce$pseudotime_bin==0] = 1
    } else {
      pseudotime_order = order(pseudotime, decreasing = FALSE)

      ncells = ncol(SingleCellExperiment::counts(pseudotime_sce))
      cells_per_bin = floor(ncells/n_bins)
      pseudotime_ordered_cells = rownames(pseudotime_sce@colData)[pseudotime_order]

      pseudotime_sce$pseudotime_bin = n_bins

      for (i in seq_len(n_bins)) {
        cells_for_bin = pseudotime_ordered_cells[(i*cells_per_bin-cells_per_bin+1) : (i*cells_per_bin)]
        pseudotime_sce[,cells_for_bin]$pseudotime_bin = i
      }
    }

    bin_ids = sort(unique(pseudotime_sce$pseudotime_bin))
    pseudobulks = as.data.frame(Matrix::rowSums(SingleCellExperiment::normcounts(pseudotime_sce)))
    for (i in bin_ids) {
      pb = as.data.frame(Matrix::rowSums(SingleCellExperiment::normcounts(subset(pseudotime_sce, , pseudotime_sce@colData[["pseudotime_bin"]] == i))))
      colnames(pb) = i
      pseudobulks = cbind(pseudobulks, pb)
    }
    pseudobulks = pseudobulks[2:length(pseudobulks)]

    return(methods::new("AtgnatData", pseudobulks = pseudobulks, bins = bin_ids))
  }
)

#' Map the best matching SC bin for a bulk sample
#'
#' @param atgnat_data The [AtgnatData] holding the bins.
#' @param bulk_id The sample id of the bulk to analyse
#' @param bulk_data The whole bulk read matrix
#'
#' @return A dataframe with one row, containing results from the mapping process.
#' * Bin: the bin that best matched the bulk sample.
#' * Corr: the spearman's rho that the test geneset had between the winning bin and the bulk.
#' * top_2_distance: the absolute difference between the best and second best mapping buckets. Higher indicates a less doubtful mapping.
#' * history: a dataframe of the correlation score for each bin.
#' @export
#'
#' @examples
map_best_bin_2 <- function(atgnat_data, bulk_id, bulk_data) {

  # TODO throw error if no genes

  sum_for_top_genes = bulk_data[atgnat_data@genes,bulk_id]

  best_cor = -1
  best_i = 0
  correlations_history = data.frame()
  for (i in atgnat_data@bins ) {
    bin_ratios = atgnat_data@pseudobulks[atgnat_data@genes,i]

    corr <- stats::cor.test(bin_ratios, sum_for_top_genes, method = 'spearman')
    if (corr$estimate > best_cor) {
      best_cor = unname(corr$estimate)
      best_i = i
    }
    correlations_history <- rbind(correlations_history, c(i, unname(corr$estimate)))
    i=i+1
  }
  colnames(correlations_history) = c('bin', 'correlation')

  top2 = utils::head(sort(correlations_history$correlation, decreasing=TRUE),n=2)
  # TODO round this to 4 decimal places not significant figures
  distance_between_top_2_corrs = signif(top2[1]-top2[2],2)

  result = data.frame(bin=best_i, correlation=best_cor, top_2_distance=distance_between_top_2_corrs, history=correlations_history)
  return(result)
}

#' Evaluate n_bins and n_genes for bin mapping
#'
#' @description Will use the n_bins and n_genes implied by the `sce` and
#' `pseudotime_bins_top_n_genes_df` parameters and return quality metrics and
#' an optional chart.
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
#' @examples
evaluate_parameters_2 <- function(atgnat_data, make_plot=FALSE, plot_columns=4) {

  bin_ids = atgnat_data@bins

  results.best_bin = c()
  results.best_corr = c()
  results.history = c()
  results.specificity = c()

  for (i in bin_ids) {
    res = map_best_bin_2(atgnat_data, i, atgnat_data@pseudobulks)
    results.best_bin = append(results.best_bin, c(res[[1,1]]))
    results.best_corr = append(results.best_corr, c(res[[1,2]]))
    results.specificity = append(results.specificity, c(res[[1,3]]))
    results.history = append(results.history, c(res[,4:5]))
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
# TODO: consider allowing this to take just a list of genes not an association test
#   result to make it a bit more reusable. If we ask it to be ordered, we can
#   just take the top n.
find_best_params_2 <- function(x, genelist, bins_count_range=c(5,10,15,20,25,30,35), gene_count_range=c(20,30,40,45,50,55,60,70,80), ...) {

  results = data.frame(gene_count=c(), bin_count=c(), worst_specificity=c(), mean_specificity=c())
  # TODO consider parallelising this
  # https://www.bioconductor.org/packages/devel/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.html#parallel-looping-vectorized-and-aggregate-operations
  # https://cran.r-project.org/web/packages/foreach/index.html
  for (bin_count in bins_count_range) {
    print(bin_count)
    atgnat_data = as.AtgnatData(x=x, n_bins=bin_count, ...)

    for (genes_count in gene_count_range) {
      atgnat_data@genes = genelist[1:genes_count]
      res = evaluate_parameters_2(atgnat_data, make_plot=FALSE)
      results = rbind(results, data.frame(bin_count=c(bin_count), gene_count=c(genes_count), worst_specificity=c(res[1]), mean_specificity=c(res[2])))
    }
  }

  return(results)

}

#' Evaluate Top Genes
#'
#' Shows plots over bins of expression of the top n genes. This is designed to help
#' identify if you have selected genes that vary over the pseudotime you have chosen
#' bins to exist over. Uses the normcounts of the SCE.
#'
#' @param atgnat_data The [AtgnatData] to get bins and expression from.
#' @param n_genes_to_plot The number of genes to plot.
#' @param plot_columns The number of columns to plot the grid with. Best as a
#' divisor of `n_genes_to_plot`.
#'
#' @export
#'
#' @examples
evaluate_top_n_genes_2 <- function(atgnat_data, n_genes_to_plot=16, plot_columns=4) {

  # TODO check genes exist

  plots = list()
  for (i in seq_len(n_genes_to_plot)) {
    plots[[i]] = PRIVATE_plot_gene_over_bins(atgnat_data@pseudobulks, atgnat_data@genes[i])
  }

  gridExtra::grid.arrange(
    top=grid::textGrob(paste(length(atgnat_data@genes), "genes"),gp=grid::gpar(fontsize=20,font=3)),
    grobs=plots,
    ncol=plot_columns
  )

}

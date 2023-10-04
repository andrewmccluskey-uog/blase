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
# TODO make this return an object, and generics to it to let it plot on an SCE (and one day a seurat object)
map_best_bin <- function(atgnat_data, bulk_id, bulk_data) {

  if (atgnat_data@genes == NA || length(atgnat_data@genes) == 0) {
    stop("No genes to map with. Please add something to the atgnat_data@genes slot.")
  }

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
  distance_between_top_2_corrs = round(top2[1]-top2[2], 4)

  result = data.frame(bin=best_i, correlation=best_cor, top_2_distance=distance_between_top_2_corrs, history=correlations_history)
  return(result)
}

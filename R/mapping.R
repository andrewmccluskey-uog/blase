#' Map the best matching SC bin for a bulk sample
#'
#' @concept mapping
#'
#' @param atgnat_data The [AtgnatData] holding the bins.
#' @param bulk_id The sample id of the bulk to analyse
#' @param bulk_data The whole bulk read matrix
#' @param calculate_pvalues Whether or not to calculate confidence values for the correlations
#'
#' @import methods
#'
#' @return A [MappingResult] object.
#' @export
#'
#' @inherit MappingResult-class examples
map_best_bin <- function(atgnat_data, bulk_id, bulk_data, calculate_pvalues=TRUE) {

  if (any(is.na(atgnat_data@genes)) || length(atgnat_data@genes) == 0) {
    stop("No genes to map with. Please add something to the atgnat_data@genes slot.")
  }

  sum_for_top_genes = bulk_data[atgnat_data@genes,bulk_id]

  best_cor = -1
  best_i = 0
  correlations_history = data.frame()
  for (i in atgnat_data@bins ) {
    bin_ratios = atgnat_data@pseudobulks[atgnat_data@genes,i]

    corr <- stats::cor.test(bin_ratios, sum_for_top_genes, method = 'spearman',exact=calculate_pvalues)
    if (corr$estimate > best_cor) {
      best_cor = unname(corr$estimate)
      best_i = i
    }
    correlations_history <- rbind(correlations_history, c(i, unname(corr$estimate), corr$p.value))
    i=i+1
  }

  colnames(correlations_history) = c('bin', 'correlation', 'pvalue')
  correlations_history$adj_pvalue = correlations_history$pvalue * nrow(correlations_history)

  top2 = utils::head(sort(correlations_history$correlation, decreasing=TRUE),n=2)
  distance_between_top_2_corrs = round(top2[1]-top2[2], 4)

  return(methods::new("MappingResult",
                      bulk_name=bulk_id,
                      best_bin=best_i,
                      best_correlation=best_cor,
                      best_pvalue=correlations_history[best_i,"pvalue"],
                      best_adj_pvalue=correlations_history[best_i,"adj_pvalue"],
                      top_2_distance=distance_between_top_2_corrs,
                      history=correlations_history)
         )

}


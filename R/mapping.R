#' Map the best matching SC bin for a bulk sample
#'
#' @concept mapping
#'
#' @param atgnat_data The [AtgnatData] holding the bins.
#' @param bulk_id The sample id of the bulk to analyse
#' @param bulk_data The whole bulk read matrix
#'
#' @import methods
#'
#' @return A [MappingResult] object.
#' @export
#'
#' @examples
#' counts_matrix <- matrix(rep(1, cells*genes), ncol=cells, nrow=genes)
#' sce <- SingleCellExperiment::SingleCellExperiment(assays=list(normcounts=counts_matrix))
#' colnames(sce) = seq_len(cells)
#'
#' atgnat_data = as.AtgnatData(sce)
#'
#' bulk_counts = matrix(rep(1, 3*genes), ncol=3, nrow=genes)
#' colnames(bulk_counts) = c("A", "B", "C")
#' result = map_best_bin(atgnat_data, "B", bulk_counts)
#' result
#' rdg
map_best_bin <- function(atgnat_data, bulk_id, bulk_data) {

  if (any(is.na(atgnat_data@genes)) || length(atgnat_data@genes) == 0) {
    stop("No genes to map with. Please add something to the atgnat_data@genes slot.")
  }

  sum_for_top_genes = bulk_data[atgnat_data@genes,bulk_id]

  best_cor = -1
  best_i = 0
  correlations_history = data.frame()
  for (i in atgnat_data@bins ) {
    bin_ratios = atgnat_data@pseudobulks[atgnat_data@genes,i]

    # TODO we probably only want to use exact=FALSE when doing the tuning calls as the lists will be identical,
    # consider reporting confidence too
    corr <- stats::cor.test(bin_ratios, sum_for_top_genes, method = 'spearman',exact=FALSE)
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

  return(methods::new("MappingResult",
                      bulk_name=bulk_id,
                      best_bin=best_i,
                      best_correlation=best_cor,
                      top_2_distance=distance_between_top_2_corrs,
                      history=correlations_history)
         )

}


#' Map the best matching SC bin for a bulk sample
#'
#' @concept mapping
#'
#' @param blase_data The [BlaseData] holding the bins.
#' @param bulk_id The sample id of the bulk to analyse.
#' @param bulk_data The whole bulk read matrix.
#' @param bootstrap_iterations The number of bootstrapping iterations to run.
#'
#' @import methods
#'
#' @return A [MappingResult] object.
#' @export
#'
#' @inherit MappingResult-class examples
map_best_bin <- function(blase_data, bulk_id, bulk_data, bootstrap_iterations=200) {

  if (any(is.na(blase_data@genes)) || length(blase_data@genes) == 0) {
    stop("No genes to map with. Please add something to the blase_data@genes slot.")
  }

  best_cor = -1
  best_i = 0
  correlations_history = data.frame()

  for (i in blase_data@bins) {

    genes_present = blase_data@genes[blase_data@genes %in% rownames(blase_data@pseudobulk_bins[[i]])]
    counts_for_top_genes = bulk_data[genes_present,as.character(bulk_id)]

    if (any(length(genes_present) != length(blase_data@genes))) {
      warn(paste('Not all genes present in bucket',i,"continuing without checking correlation for these genes.\n"))
    }
    
    
    #print(dim(blase_data@pseudobulk_bins[[i]]))
    if (ncol(blase_data@pseudobulk_bins[[i]]) <= 1) {
      stop(paste0('Not enough cells in bin ', i,' to map against, please reduce number of bins (currently ', length(blase_data@pseudobulk_bins),') or split by cells'))
    }

    bin_ratios = blase_data@pseudobulk_bins[[i]][genes_present,]

    all_info_correlation = stats::cor.test(unname(Matrix::rowMeans(bin_ratios)), counts_for_top_genes, method='spearman', exact=FALSE)
    corr_estimate = unname(all_info_correlation$estimate)

    lower_bound = 0
    upper_bound = 1

    if (bootstrap_iterations > 0) {
      correlations = c()
      for (j in seq_len(bootstrap_iterations)) {
        pseudobulk_sample_names = sample(colnames(bin_ratios), ncol(bin_ratios), replace = TRUE)
        pseudobulk_sample_means = unname(Matrix::rowMeans(bin_ratios[,pseudobulk_sample_names]))

        corr <- stats::cor.test(pseudobulk_sample_means, counts_for_top_genes, method = 'spearman',exact=FALSE)
        correlations = c(correlations, unname(corr$estimate))
      }

      middle_90 = correlations[{q<-rank(correlations)/length(correlations);q<0.05 | q>=0.95}]
      lower_bound = min(middle_90)
      upper_bound = max(middle_90)
    }

    if (corr_estimate > best_cor) {
      best_cor = corr_estimate
      best_i = i
    }
    correlations_history <- rbind(correlations_history, c(i, corr_estimate, lower_bound, upper_bound))
    i=i+1
  }

  colnames(correlations_history) = c('bin', 'correlation', 'lower_bound', 'upper_bound')

  top2 = utils::head(sort(correlations_history$correlation, decreasing=TRUE),n=2)
  distance_between_top_2_corrs = round(top2[1]-top2[2], 4)

  confident_mapping = FALSE
  non_top_mapping_best_upper_bound = max(
    correlations_history[correlations_history$bin!=best_i,]$upper_bound
  )
  top_mapping_lower_bound = correlations_history[correlations_history$bin==best_i,]$lower_bound

  ## TODO surface to the user why the mapping wasn't confident, and allow user override of the cutoff
  confident_mapping = non_top_mapping_best_upper_bound < top_mapping_lower_bound && top_mapping_lower_bound > 0

  return(methods::new("MappingResult",
                      bulk_name=bulk_id,
                      best_bin=best_i,
                      best_correlation=best_cor,
                      top_2_distance=distance_between_top_2_corrs,
                      confident_mapping=confident_mapping,
                      history=correlations_history,
                      bootstrap_iterations=bootstrap_iterations)
         )

}


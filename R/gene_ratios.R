#' Get Top Genes From An AssociationTestResult
#'
#' Pulls the genes with the highest wald statistic from an association test result, with a p value cutoff.
#' @param association_test_results The association test results data frame to take the genes from.
#' @param n_genes The number of genes to return. Defaults to 40.
#' @param lineage The Lineage to use. The Defaults to NA, which assumes the test was run with Lineages=False.
#' @param p_cutoff The P value cutoff to use. Defaults to less than 0.05.
#'
#' @return A vector of the names of the genes that best describe a lineage's trajectory.
#' @export
#'
#' @examples
get_top_n_genes = function(association_test_results, n_genes=40, lineage=NA, p_cutoff=0.05) {

  pvalue_slot_for_lineage = "pvalue"
  wald_slot_for_lineage = "waldStat"

  if (!is.na(lineage)) {
    pvalue_slot_for_lineage = paste0("pvalue_", lineage)
    wald_slot_for_lineage = paste0("waldStat_", lineage)
  }

  # P Cutoff
  asso_results_copy = association_test_results[association_test_results[,pvalue_slot_for_lineage] < p_cutoff,]
  # Remove NAs
  asso_results_copy = asso_results_copy[!is.na(rownames(asso_results_copy)), ]
  # Sort by wald stat
  asso_results_copy = asso_results_copy[order(-asso_results_copy[,wald_slot_for_lineage]), ]

  topPdtGenesNames <- rownames(asso_results_copy)[1:n_genes]
  return(topPdtGenesNames)

}


#' Get a dataframe of gene expression for a gene list in all bins
#'
#' @param pseudotime_sce The SingleCellExperiment to get the counts from.
#' @param gene_list The list of genes to summarise.
#'
#' @return A dataframe of total gene expression for the selected genes by pseudotime bin.
#' @export
#'
#' @examples
get_pseudotime_gene_ratios = function(pseudotime_sce, gene_list) {

  bin_ids = as.numeric(names(table(pseudotime_sce@colData["pseudotime_bin"])))

  pseudotime_bins_ratios_of_top_n_genes_df = data.frame()
  for (i in bin_ids) {
    sce_of_bin_n = subset(pseudotime_sce, , pseudotime_bin==i)

    ## Get the mean of the expression of each gene at this pseudotime bin
    mean_for_bin_top_genes = rowSums(SingleCellExperiment::logcounts(sce_of_bin_n[gene_list,]))

    ## get expression of genes for cells in each bin
    proportion_for_bin_top_genes = rowSums(SingleCellExperiment::logcounts(sce_of_bin_n[gene_list,]))
    pseudotime_bins_ratios_of_top_n_genes_df = rbind(pseudotime_bins_ratios_of_top_n_genes_df, proportion_for_bin_top_genes)
  }
  colnames(pseudotime_bins_ratios_of_top_n_genes_df) = gene_list
  return(pseudotime_bins_ratios_of_top_n_genes_df)

}

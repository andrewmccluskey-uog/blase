#' Get Top Genes From An AssociationTestResult
#'
#' Pulls the genes with the highest wald statistic from an association test result, with a p value cutoff.
#'
#' @concept util
#'
#' @param association_test_results The association test results data frame to take the genes from.
#' @param n_genes The number of genes to return. Defaults to 40.
#' @param lineage The Lineage to use. The Defaults to NA, which assumes the test was run with Lineages=False.
#' @param p_cutoff The P value cutoff to use. Defaults to less than 0.05.
#'
#' @return A vector of the names of the genes that best describe a lineage's trajectory.
#' @export
#'
#' @examples
#' assoRes = data.frame(
#'   row.names=c("A", "B", "C", "D"),
#'   waldStat=c(25, 50, 100, 10),
#'   pvalue=c(0.01, 0.5, 0.005, 0.13))
#' get_top_n_genes(assoRes, n_genes=2)
get_top_n_genes <- function(association_test_results, n_genes=40, lineage=NA, p_cutoff=0.05) {

  pvalue_slot_for_lineage <- "pvalue"
  wald_slot_for_lineage <- "waldStat"

  if (!is.na(lineage)) {
    pvalue_slot_for_lineage <- paste0("pvalue_", lineage)
    wald_slot_for_lineage <- paste0("waldStat_", lineage)
  }

  # P Cutoff
  asso_results_copy <- association_test_results[association_test_results[,pvalue_slot_for_lineage] < p_cutoff,]
  # Remove NAs
  asso_results_copy <- asso_results_copy[!is.na(rownames(asso_results_copy)), ]
  # Sort by wald stat
  asso_results_copy <- asso_results_copy[order(-asso_results_copy[,wald_slot_for_lineage]), ]

  topPdtGenesNames <- rownames(asso_results_copy)[seq_len(n_genes)]
  topPdtGenesNames <- topPdtGenesNames[!is.na(topPdtGenesNames)]
  return(topPdtGenesNames)

}

#' Get a pseudobulk of bins with at least 2 replicates
#'
#' This function will try to create a pseudobulked count matrix for the bins.
#' When a replicate has too few cells, it is discounted. If only one exists,
#' then we sample from it twice to create the pseudobulks.
#'
#' @concept util
#'
#' @param pseudotime_sce The Single Cell Experiment object to get the bins from
#' @param min_cells_for_bulk The minimum cells to look for per replicate and bin.
#' @param replicate_slot The slot in the Single Cell Experiment that contains replicate information
#'
#' @return A dataframe pseudobulk counts matrix.
#' @export
#'
#' @examples
#' library(SingleCellExperiment, quietly=TRUE)
#' library(blase)
#' counts <- matrix(rpois(1000, lambda = 10), ncol=100, nrow=10)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(normcounts = counts, counts = counts/2)
#' )
#' sce$pseudotime = seq_len(100)
#' colnames(sce) = seq_len(100)
#' rownames(sce) = as.character(seq_len(10))
#' sce = assign_pseudotime_bins(sce, n_bins=5,
#'   pseudotime_slot="pseudotime", split_by="cells")
#' sce$replicate=rep(c(1,2), 50)
#' result = get_bins_as_bulk(sce, min_cells_for_bulk=1, replicate_slot="replicate")
#' result
get_bins_as_bulk <- function(pseudotime_sce, min_cells_for_bulk=50, replicate_slot="replicate") {

  # TODO Generalize so we don't rely on subsetting with subset on a magic field (i.e. replicate, which is even listed as a param but isn't really)

  output <- data.frame()
  for (bin_id in seq_len(max(pseudotime_sce$pseudotime_bin))) {

    bin_specific_sce <- subset(pseudotime_sce, , pseudotime_sce@colData[["pseudotime_bin"]] == bin_id)

    counts <- table(bin_specific_sce@colData[[replicate_slot]])

    replicates_with_more_than_minimum <- rownames(as.data.frame(counts[counts>min_cells_for_bulk]))

    if (length(replicates_with_more_than_minimum) >= 2) {

      pseudobulks <- data.frame()
      for (rep_id in replicates_with_more_than_minimum) {

        bin_specific_rep_specific_sce_pseudobulk <- as.data.frame(rowSums(SingleCellExperiment::counts(subset(bin_specific_sce, , bin_specific_sce@colData[["replicate"]] == rep_id))))

        if (ncol(pseudobulks) == 0) {
          pseudobulks <- bin_specific_rep_specific_sce_pseudobulk
          colnames(pseudobulks) <- c(paste0("bin_",bin_id, "_rep_", rep_id))
        } else {
          pseudobulks <- merge(x=pseudobulks, y=bin_specific_rep_specific_sce_pseudobulk, by="row.names")
          rownames(pseudobulks) <- pseudobulks$Row.names
          pseudobulks <- pseudobulks[,c(-1)]
          colnames(pseudobulks) <- append(colnames(pseudobulks)[seq_len(ncol(pseudobulks)-1)], paste0("bin_",bin_id, "_rep_", rep_id))
        }

      }


    } else if (length(replicates_with_more_than_minimum) == 1) {

      bin_specific_pseudobulk <- subset(bin_specific_sce, , bin_specific_sce@colData[["replicate"]] == replicates_with_more_than_minimum[1])

      cell_count_to_sample <- ceiling((ncol(SingleCellExperiment::counts(bin_specific_pseudobulk)) * 0.75))
      counts_for_bulk <- SingleCellExperiment::counts(bin_specific_pseudobulk)
      pseudobulks <- merge(
        x=as.data.frame(rowSums(counts_for_bulk[,sample(ncol(counts_for_bulk), size=cell_count_to_sample)])),
        y=as.data.frame(rowSums(counts_for_bulk[,sample(ncol(counts_for_bulk), size=cell_count_to_sample)])),
        by="row.names")
      rownames(pseudobulks) <- pseudobulks$Row.names
      pseudobulks <- pseudobulks[,c(-1)]
      colnames(pseudobulks) <- c(paste0("bin_",bin_id, "_rep_", replicates_with_more_than_minimum[1], "_1"), paste0("bin_",bin_id, "_rep_", replicates_with_more_than_minimum[1], "_2"))
    } else {
      message("Couldn't create pseudobulks due to too few cells in every replicate for bin ", bin_id)
      next()
    }

    if (ncol(output) == 0) {
      output <- pseudobulks
    } else {
      output <- merge(output, pseudobulks, by="row.names")
      rownames(output) <- output$Row.names
      output <- output[,c(-1)]
    }

  }

  return(output)

}

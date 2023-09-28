#' Get a pseudobulk of bins with at least 2 replicates
#'
#' This function will try to create a pseudobulked count matrix for the bins.
#' When a replicate has too few cells, it is discounted. If only one exists,
#' then we sample from it twice to create the pseudobulks.
#'
#' @param pseudotime_sce The Single Cell Experiment object to get the bins from
#' @param min_cells_for_bulk The minimum cells to look for per replicate and bin.
#' @param replicate_slot The slot in the Single Cell Experiment that contains replicate information
#'
#' @return A dataframe pseudobulk counts matrix.
#' @export
#'
#' @examples
# library(SingleCellExperiment, quietly=TRUE)
# library(atgnat)
# counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
# sce <- SingleCellExperiment::SingleCellExperiment(counts)
# sce$pseudotime = seq_len(10)
# sce = create_pseudotime_bins(sce, 5, pseudotime_slot="pseudotime")
# sce$replicate=rep(c(1,2), 5)
# result = get_bins_as_bulk(sce, min_cells_for_bulk=1, replicate_slot="replicate")
# result
get_bins_as_bulk <- function(pseudotime_sce, min_cells_for_bulk=50, replicate_slot="replicate") {

  # TODO Generalize so we don't rely on subsetting with subset on a magic field (i.e. replicate, which is even listed as a param but isn't really)

  output = data.frame()

  for (bin_id in seq_len(max(pseudotime_sce$pseudotime_bin))) {
    bin_specific_sce = subset(pseudotime_sce, , pseudotime_sce@colData[["pseudotime_bin"]] == bin_id)

    counts = table(bin_specific_sce@colData[replicate_slot])
    replicates_with_more_than_minimum = rownames(as.data.frame(counts[counts>min_cells_for_bulk]))

    if (length(replicates_with_more_than_minimum) >= 2) {

      pseudobulks = data.frame()
      for (rep_id in replicates_with_more_than_minimum) {

        bin_specific_rep_specific_sce_pseudobulk = as.data.frame(rowSums(SingleCellExperiment::counts(subset(bin_specific_sce, , bin_specific_sce@colData[["replicate"]] == rep_id))))

        if (ncol(pseudobulks) == 0) {
          pseudobulks = bin_specific_rep_specific_sce_pseudobulk
          colnames(pseudobulks) = c(paste0("bin_",bin_id, "_rep_", rep_id))
        } else {
          pseudobulks = merge(x=pseudobulks, y=bin_specific_rep_specific_sce_pseudobulk, by="row.names")
          rownames(pseudobulks) = pseudobulks$Row.names
          pseudobulks = pseudobulks[,c(-1)]
          colnames(pseudobulks) = append(colnames(pseudobulks)[1:ncol(pseudobulks)-1], paste0("bin_",bin_id, "_rep_", rep_id))
        }

      }


    } else if (length(replicates_with_more_than_minimum) == 1) {

      bin_specific_pseudobulk = subset(bin_specific_sce, , bin_specific_sce@colData[["replicate"]] == replicates_with_more_than_minimum[1])

      cell_count_to_sample = ceiling((ncol(SingleCellExperiment::counts(bin_specific_pseudobulk)) * 0.75))
      counts_for_bulk = SingleCellExperiment::counts(bin_specific_pseudobulk)
      pseudobulks = merge(
        x=as.data.frame(rowSums(counts_for_bulk[,sample(ncol(counts_for_bulk), size=cell_count_to_sample)])),
        y=as.data.frame(rowSums(counts_for_bulk[,sample(ncol(counts_for_bulk), size=cell_count_to_sample)])),
        by="row.names")
      rownames(pseudobulks) = pseudobulks$Row.names
      pseudobulks = pseudobulks[,c(-1)]
      colnames(pseudobulks) = c(paste0("bin_",bin_id, "_rep_", replicates_with_more_than_minimum[1], "_1"), paste0("bin_",bin_id, "_rep_", replicates_with_more_than_minimum[1], "_2"))
    } else {
      print(paste("Couldn't create pseudobulks due to too few cells in every replicate ", bin_id))
      next()
    }

    if (ncol(output) == 0) {
      output = pseudobulks
    } else {
      output = merge(output, pseudobulks, by="row.names")
      rownames(output) = output$Row.names
      output = output[,c(-1)]
    }

  }

  return(output)

}

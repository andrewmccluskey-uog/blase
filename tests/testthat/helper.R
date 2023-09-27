generate_test_sce = function(cells=100, genes=50) {
  counts_matrix <- matrix(rpois(cells*genes, lambda = 10), ncol=cells, nrow=genes)
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts_matrix))
  colnames(sce) = seq_len(cells)
  return(sce)
}

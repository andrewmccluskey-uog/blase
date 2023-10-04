generate_test_sce = function(cells=100, genes=50) {
  counts_matrix <- matrix(rep(1, cells*genes), ncol=cells, nrow=genes)
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(normcounts=counts_matrix))
  colnames(sce) = seq_len(cells)
  return(sce)
}

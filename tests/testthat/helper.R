generate_test_sce = function(cells=100, genes=50) {
  counts_matrix <- matrix(rep(1, cells*genes), ncol=cells, nrow=genes)
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts_matrix*3, normcounts=counts_matrix))
  colnames(sce) = seq_len(cells)
  rownames(sce) = seq_len(genes)
  return(sce)
}

generate_test_atgnat_data = function(cells=100, genes=50) {
  sce = generate_test_sce(cells, genes)
  atgnat_data = as.AtgnatData(sce)
  atgnat_data@genes = seq_len(genes)
  return(atgnat_data)
}

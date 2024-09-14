generate_test_sce <- function(cells = 100, genes = 50) {
    counts_matrix <- matrix(rep(1, cells * genes), ncol = cells, nrow = genes)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_matrix * 3, normcounts = counts_matrix, logcounts = counts_matrix))
    sce$pseudotime <- (1:cells) / cells
    colnames(sce) <- seq_len(cells)
    rownames(sce) <- seq_len(genes)
    return(sce)
}

generate_test_seurat <- function(cells = 100, genes = 50) {
    sce <- generate_test_sce(cells, genes)
    return(Seurat::as.Seurat(sce))
}

generate_test_blase_data <- function(cells = 100, genes = 50) {
    sce <- generate_test_sce(cells, genes)
    blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 5)
    genes(blase_data) <- as.character(seq_len(genes))
    return(blase_data)
}

generate_test_sce <- function(cells = 100, genes = 50) {
    counts_matrix <- matrix(rep(1, cells * genes), ncol = cells, nrow = genes)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_matrix * 3, normcounts = counts_matrix, logcounts = counts_matrix))
    sce$pseudotime <- (1:cells) / cells
    colnames(sce) <- paste0("C", seq_len(cells))
    rownames(sce) <- paste0("G", seq_len(genes))
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

# Underscores below because some of these are function names exposed by BLASE
generate_test_mapping_result <- function(
    name_ = "Test",
    best_bin_ = 1,
    best_correlation_ = 0.05,
    top_2_distance_ = 0.001,
    confident_mapping_ = FALSE,
    bootstrap_iterations_ = 200,
    history_ = NULL) {
    if (is.null(history_)) {
        history_ <- data.frame(
            bin = c(1, 2),
            correlation = c(0.5, 0.2),
            lower_bound = c(0.3, 0.1),
            upper_bound = c(0.6, 0.3)
        )
    }
    return(methods::new("MappingResult",
        bulk_name = name_,
        best_bin = best_bin_,
        best_correlation = best_correlation_,
        top_2_distance = top_2_distance_,
        confident_mapping = confident_mapping_,
        bootstrap_iterations = bootstrap_iterations_,
        history = history_
    ))
}

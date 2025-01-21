## evaluate_parameters
test_that("evaluate_parameters() generates tuple of minimum and mean specifcity, pct confident mappings for parameters", {
    cells <- 100
    genes <- 20
    counts_matrix <- matrix(seq_len(cells * genes), ncol = cells, nrow = genes)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_matrix * 3, normcounts = counts_matrix, logcounts = counts_matrix / 2))
    sce$pseudotime <- (1:cells) / cells
    colnames(sce) <- paste0("C", seq_len(cells))
    rownames(sce) <- paste0("G", seq_len(genes))

    blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 5)
    blase_data@genes <- paste0("G", seq_len(genes))

    result <- evaluate_parameters(blase_data)

    expect_equal(result, c(0, 0, 0))
})

## find_best_params
test_that("find_best_params() warns about too few genes", {
    sce <- generate_test_sce(200, 300)

    tmp1 <- function() find_best_params(sce, genelist = seq_len(50))
    expect_error(tmp1(), "Not enough genes provided to meet tuning requests. Provided=50 wanted=80", fixed = TRUE)
})

test_that("find_best_params() runs without error", {
  sce <- generate_test_sce(200, 300)
  sce$pseudotime <- (1:200) / 200
  colnames(sce) <- paste0("C", seq_len(200))
  rownames(sce) <- paste0("G", seq_len(300))

  results = find_best_params(
    sce,
    genelist = paste0("G",seq_len(40)),
    bins_count_range = 2:3,
    gene_count_range = c(5,10),
    pseudotime_slot = "pseudotime")

  expected = data.frame(data.frame(
    column_label=c("1", "2", "1", "2"),
    bin_count=c(2, 2, 3, 3),
    gene_count=c(5, 10, 5, 10),
    min_convexity=as.double(c(NA, NA, NA, NA)),
    mean_convexity=as.double(c(NA, NA, NA, NA)),
    confident_mapping_pct=as.double(c(NA, NA, NA,NA))
  ))
  rownames(expected) = seq_len(4)

  expect_equal(results, expected)
})

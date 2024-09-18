## evaluate_parameters
test_that("evaluate_parameters() generates tuple of worst and mean specifcity for parameters", {
  cells <- 100
  genes <- 20
  counts_matrix <- matrix(seq_len(cells * genes), ncol = cells, nrow = genes)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_matrix * 3, normcounts = counts_matrix, logcounts = counts_matrix / 2))
  sce$pseudotime <- (1:cells) / cells
  colnames(sce) <- seq_len(cells)
  rownames(sce) <- seq_len(genes)

  blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 5)
  blase_data@genes <- as.character(seq_len(genes))

  result <- evaluate_parameters(blase_data)

  expect_equal(result, c(0, 0))
})

## find_best_params
test_that("find_best_params() warns about too few genes", {
  sce <- generate_test_sce(200, 300)

  tmp1 <- function() find_best_params(sce, genelist = seq_len(50))
  expect_error(tmp1(), "Not enough genes provided to meet tuning requests. Provided=50 wanted=80", fixed = TRUE)
})

## TODO testing ggplots = https://stackoverflow.com/questions/31038709/how-to-write-a-test-for-a-ggplot-plot

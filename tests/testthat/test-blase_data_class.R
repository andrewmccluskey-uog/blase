## Test genes accessors
test_that("genes can be set using accessor", {
  cells <- 15
  genes <- 20

  sce <- generate_test_sce(cells, genes)
  blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 5)
  genes(blase_data) <- as.character(seq_len(genes))

  expect_equal(blase_data@genes, as.character(seq_len(genes)))
})

test_that("genes can be gotten using accessor", {
  cells <- 15
  genes <- 20

  sce <- generate_test_sce(cells, genes)
  blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 5)
  blase_data@genes <- as.character(seq_len(genes))

  expect_equal(genes(blase_data), as.character(seq_len(genes)))
})

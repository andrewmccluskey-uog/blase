test_that("generates correct number of bins", {
  sce = generate_test_sce(cells=150, genes=100)
  sce$pseudotime = sqrt(seq_len(150))

  atgnat_data = as.AtgnatData(sce, pseudotime_slot="pseudotime", n_bins=5, split_by="pseudotime_range")
  expect_equal(length(atgnat_data@bins), 5)
  expect_equal(ncol(atgnat_data@pseudobulks), 5)
})


test_that("generates pseudotime bins by pseudotime_range", {
  n_cells = 150
  sce = generate_test_sce(cells=n_cells, genes=100)
  sce$pseudotime = sqrt(seq_len(n_cells))

  atgnat_data = as.AtgnatData(sce, pseudotime_slot="pseudotime", n_bins=6, split_by="pseudotime_range")

  # In tests, gene expression per cell == 1
  cells_per_bin = colMeans(atgnat_data@pseudobulks)
  expect_equal(n_cells, sum(colMeans(atgnat_data@pseudobulks)))
  expected_output =  c("1"=4, "2"=14, "3"=24, "4"=33, "5"=42, "6"=33)
  expect_equal(cells_per_bin, expected_output)
})

test_that("generates pseudotime bins by cells", {
  n_cells = 150
  sce = generate_test_sce(cells=n_cells, genes=100)
  sce$pseudotime = sqrt(seq_len(n_cells))

  atgnat_data = as.AtgnatData(sce, pseudotime_slot="pseudotime", n_bins=6, split_by="cells")

  cells_per_bin = colMeans(atgnat_data@pseudobulks) # In tests, gene expression per cell == 1
  expect_equal(n_cells, sum(colMeans(atgnat_data@pseudobulks)))
  expected_output = c("1"=25, "2"=25, "3"=25, "4"=25, "5"=25, "6"=25)
  expect_equal(cells_per_bin, expected_output)
})

test_that("throws error for invalid split_by parameter", {
  sce = generate_test_sce(cells=150, genes=100)
  sce$pseudotime = sqrt(seq_len(150))

  tmp1 <- function() atgnat_data = as.AtgnatData(sce, pseudotime_slot="pseudotime", n_bins=10, split_by="this_won't_work")
  expect_error(tmp1(), "split_by must be 'pseudotime_range' or 'cells'", fixed=TRUE)
})



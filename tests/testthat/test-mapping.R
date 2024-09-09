## WIP test map_best_bin

test_that("throws error if gene list is null", {
  
  sce = generate_test_sce()
  blase_data = as.BlaseData(sce, pseudotime_slot="pseudotime", n_bins=5)
  
  tmp1 <- function() map_best_bin(blase_data, "1", data.frame(), bootstrap_iterations=1)
  expect_error(tmp1(), "No genes to map with. Please add something to the blase_data@genes slot.", fixed=TRUE)
  
})

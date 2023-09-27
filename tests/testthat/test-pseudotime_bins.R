test_that("generates pseudotime bins by pseudotime_range", {
  sce = generate_test_sce(cells=500, genes=100)
  sce$pseudotime = sqrt(seq_len(500))

  sce = create_pseudotime_bins(sce, 10, "pseudotime", split_by="pseudotime_range")
  result = as.data.frame(table(sce$pseudotime_bin, useNA="always"))

  expected_output =  data.frame(Var1=c(1,2,3,4,5,6,7,8,9,10, NA), Freq=c(5,16,26,37,48,58,69,79,90,72, 0))
  expected_output$Var1 = as.factor(expected_output$Var1)

  expect_equal(result, expected_output)
})

test_that("generates pseudotime bins by cells", {
  sce = generate_test_sce(cells=500, genes=100)
  sce$pseudotime = sqrt(seq_len(500))

  sce = create_pseudotime_bins(sce, 10, "pseudotime", split_by="cells")

  result = as.data.frame(table(sce$pseudotime_bin, useNA="always"))

  expected_output =  data.frame(Var1=c(1,2,3,4,5,6,7,8,9,10, NA), Freq=c(50,50,50,50,50,50,50,50,50,50, 0))
  expected_output$Var1 = as.factor(expected_output$Var1)

  expect_equal(result, expected_output)
})

test_that("throws error for invalid split_by parameter", {
  sce = generate_test_sce(cells=500, genes=100)
  sce$pseudotime = sqrt(seq_len(500))

  tmp1 <- function() create_pseudotime_bins(sce, 10, "pseudotime", split_by="blah")
  expect_error(tmp1(), "split_by must be 'pseudotime_range' or 'cells'", fixed=TRUE)
})



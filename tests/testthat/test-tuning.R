test_that("find_best_params() warns about too few genes", {
  sce = generate_test_sce(200, 300)

  tmp1 <- function() find_best_params(sce, genelist = seq_len(50))
  expect_error(tmp1(), "Not enough genes provided to meet tuning requests. Provided=50 wanted=80", fixed=TRUE)
})

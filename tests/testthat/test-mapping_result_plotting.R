test_that("plot_mapping_result_heatmap throws error if list has no MappingResults (char)", {
  tmp1 <- function() plot_mapping_result_heatmap(list("abc"))
  expect_error(tmp1(), "You must provide a list of MappingResult objects only.", fixed = TRUE)
})

test_that("plot_mapping_result_heatmap throws error if list has no MappingResults (object)", {
  tmp1 <- function() plot_mapping_result_heatmap(list(SingleCellExperiment::SingleCellExperiment()))
  expect_error(tmp1(), "You must provide a list of MappingResult objects only.", fixed = TRUE)
})

test_that("plot_mapping_result_heatmap throws error if list has a mix of MappingResults and other types", {
  mapping_result <- methods::new("MappingResult",
    bulk_name = "Test",
    best_bin = 1,
    best_correlation = 0.05,
    top_2_distance = 0.001,
    history = data.frame()
  )

  tmp1 <- function() plot_mapping_result_heatmap(list(mapping_result, 5))
  expect_error(tmp1(), "You must provide a list of MappingResult objects only.", fixed = TRUE)
})

test_that("plot_mapping_result_heatmap runs for list of MappingResults only", {
  mapping_result <- methods::new("MappingResult",
    bulk_name = "Test",
    best_bin = 1,
    best_correlation = 0.05,
    top_2_distance = 0.001,
    confident_mapping = TRUE,
    history = data.frame(bin = c(1), correlation = c(0.5), lower_bound = c(0.3), upper_bound = c(0.6)),
    bootstrap_iterations = 200
  )

  tmp1 <- function() plot_mapping_result_heatmap(list(mapping_result))
  expect_no_error(tmp1())
})

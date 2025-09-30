# Tests for plot_mapping_result
test_that("plot_mapping_result runs without error", {
  sce <- generate_test_sce()
  SingleCellExperiment::reducedDim(sce, "UMAP") <- uwot::umap(t(SingleCellExperiment::logcounts(sce)), n_neighbors = 2)
  sce <- assign_pseudotime_bins(sce, n_bins=3, pseudotime_slot="pseudotime")
  sce$cluster = round(runif(ncol(sce), 1, 3))
  mapping_result <- generate_test_mapping_result()

  tmp1 <- function() print(plot_mapping_result(sce, mapping_result, "cluster"))
  expect_no_error(tmp1())
})

# Tests for plot_mapping_result_heatmap
test_that("plot_mapping_result_heatmap throws error if list has no MappingResults (char)", {
    tmp1 <- function() plot_mapping_result_heatmap(list("abc"))
    expect_error(tmp1(), "You must provide a list of MappingResult objects only.", fixed = TRUE)
})

test_that("plot_mapping_result_heatmap throws error if list has no MappingResults (object)", {
    tmp1 <- function() plot_mapping_result_heatmap(list(SingleCellExperiment::SingleCellExperiment()))
    expect_error(tmp1(), "You must provide a list of MappingResult objects only.", fixed = TRUE)
})

test_that("plot_mapping_result_heatmap throws error if list has a mix of MappingResults and other types", {
    mapping_result <- generate_test_mapping_result()

    tmp1 <- function() plot_mapping_result_heatmap(list(mapping_result, 5))
    expect_error(tmp1(), "You must provide a list of MappingResult objects only.", fixed = TRUE)
})

test_that("plot_mapping_result_heatmap runs for list of MappingResults only", {
    mapping_result <- generate_test_mapping_result()

    tmp1 <- function() print(plot_mapping_result_heatmap(list(mapping_result)))
    expect_no_error(tmp1())
})

test_that("plot_mapping_result_heatmap annotates best match boxes when annotate true", {
  mapping_result <- generate_test_mapping_result()

  p = plot_mapping_result_heatmap(list(mapping_result))
  expect_s3_class(p, "ggplot")
  expect_equal(p$layers$geom_tile$aes_params$linewidth, 0.8)
  expect_equal(p$scales$scales[[2]]$breaks, c(FALSE, TRUE))

  built_p = ggplot2::ggplot_build(p)
  expect_equal(built_p$data[[1]]$colour, c("black", "transparent"))

})

test_that("plot_mapping_result_heatmap runs for list of MappingResults only", {
  mapping_result <- generate_test_mapping_result()

  p = plot_mapping_result_heatmap(list(mapping_result), annotate_confidence = FALSE)
  expect_s3_class(p, "ggplot")
  expect_equal(p$scales$scales[[2]]$breaks, c(FALSE, TRUE))

  built_p = ggplot2::ggplot_build(p)
  expect_equal(built_p$data[[1]]$colour, c("transparent", "transparent"))

})

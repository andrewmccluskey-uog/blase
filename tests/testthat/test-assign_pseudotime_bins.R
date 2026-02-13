# assign_pseudotime_bins (dataframe)
test_that("(pdt_range) throws error as pseudotime bins cannot be applied to bulks", {
    bulk_counts <- matrix(seq_len(15) * 10, ncol = 3, nrow = 5)
    colnames(bulk_counts) <- c("A", "B", "C")
    rownames(bulk_counts) <- as.character(seq_len(5))
    bulk_counts <- as.data.frame(bulk_counts)

    tmp1 <- function() assign_pseudotime_bins(bulk_counts)
    expect_error(tmp1(), "Can't update bulk data, using each sample as bins.", fixed = TRUE)
})

# assign_pseudotime_bins (seurat)
test_that("(pdt_range) adds pseudotime bins to a seurat object", {
    seurat <- generate_test_seurat()
    seurat <- assign_pseudotime_bins(seurat, "pseudotime_range", 2, pseudotime_slot = "pseudotime")
    expect_equal(unname(seurat$pseudotime_bin), c(rep(1, 50), rep(2, 50)))
})

# assign_pseudotime_bins (sce)
test_that("(pdt_range) adds pseudotime bins to a sce object", {
    sce <- generate_test_sce()
    sce <- assign_pseudotime_bins(sce, "pseudotime_range", 2, pseudotime_slot = "pseudotime")
    expect_equal(unname(sce$pseudotime_bin), c(rep(1, 50), rep(2, 50)))
})

test_that("(pdt_range) correctly places pseudotime 0 cell in bin 1", {
    sce <- generate_test_sce()
    sce$pseudotime[1] <- 0
    sce <- assign_pseudotime_bins(sce, "pseudotime_range", 2, pseudotime_slot = "pseudotime")
    expect_equal(unname(sce$pseudotime_bin), c(rep(1, 50), rep(2, 50)))
})

test_that("(pdt_range) adds correct number of bins when max pseudotime < 1", {
    sce <- generate_test_sce()
    sce$pseudotime <- seq_len(ncol(sce)) / (ncol(sce) * 2)
    sce$pseudotime[1] <- 0
    sce <- assign_pseudotime_bins(sce, "pseudotime_range", 10, pseudotime_slot = "pseudotime")
    expect_equal(length(unique(sce$pseudotime_bin)), 10)
})

test_that("(pdt_range) throws error if pseudotime starts above 0", {
  sce <- generate_test_sce()
  sce$pseudotime <- seq_len(ncol(sce)) / (ncol(sce) * 2)
  sce$pseudotime[1] <- 1

  tmp1 <- function() assign_pseudotime_bins(sce, "pseudotime_range", 10, pseudotime_slot = "pseudotime")
  expect_error(tmp1(), "pseudotime must start at 0", fixed = TRUE)
})

test_that("(pdt_range) throws error if pseudotime slot not present", {
  sce <- generate_test_sce()
  sce$pseudotime <- seq_len(ncol(sce)) / (ncol(sce) * 2)
  sce$pseudotime[1] <- 0

  tmp1 <- function() assign_pseudotime_bins(sce, "pseudotime_range", 10, pseudotime_slot = "pseudotime2")
  expect_error(tmp1(), "Pseudotime slot 'pseudotime2' does not exist", fixed = TRUE)
})

test_that("(cells) throws error if not all pseudotime slots present", {
  sce <- generate_test_sce()
  sce$pseudotime <- seq_len(ncol(sce)) / (ncol(sce) * 2)
  sce$pseudotime[1] <- 0

  tmp1 <- function() assign_pseudotime_bins(sce, "cells", 10, pseudotime_slot = c("pseudotime", "something else"))
  expect_error(tmp1(), "Pseudotime slot 'something else' does not exist", fixed = TRUE)
})

test_that("(pdt_range) throws error if nbins is 0", {
  sce <- generate_test_sce()
  sce$pseudotime <- seq_len(ncol(sce)) / (ncol(sce) * 2)
  sce$pseudotime[1] <- 0

  tmp1 <- function() assign_pseudotime_bins(sce, "pseudotime_range", 0, pseudotime_slot = "pseudotime")
  expect_error(tmp1(), "Number of bins must be greater than 0", fixed = TRUE)
})

test_that("(pdt_range) throws error if nbins is more than number of cells", {
  sce <- generate_test_sce()
  sce$pseudotime <- seq_len(ncol(sce)) / (ncol(sce) * 2)
  sce$pseudotime[1] <- 0

  tmp1 <- function() assign_pseudotime_bins(sce, "pseudotime_range", ncol(sce)+10, pseudotime_slot = "pseudotime")
  expect_error(tmp1(), "Number of bins must be less than cells in object", fixed = TRUE)
})

test_that("(pdt_range) throws error if multiple lineages", {
  sce <- generate_test_sce()
  sce$pseudotime <- seq_len(ncol(sce)) / (ncol(sce) * 2)
  sce$pseudotime[1] <- 0

  tmp1 <- function() assign_pseudotime_bins(sce, "pseudotime_range", ncol(sce)+10, pseudotime_slot = c("pseudotime"))
  expect_error(tmp1(), "Number of bins must be less than cells in object", fixed = TRUE)
})

test_that("(cell) adds pseudotime bins to a sce object", {
    sce <- generate_test_sce()
    sce <- assign_pseudotime_bins(sce, "cells", 2, pseudotime_slot = "pseudotime")
    expect_equal(unname(sce$pseudotime_bin), c(rep(1, 50), rep(2, 50)))
})

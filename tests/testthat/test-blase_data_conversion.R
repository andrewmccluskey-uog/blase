# These tests also cover assign_pseudotime_bins
test_that("generates correct number of bins", {
    sce <- generate_test_sce(cells = 150, genes = 100)
    sce$pseudotime <- sqrt(seq_len(150))

    blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 5, split_by = "pseudotime_range")
    expect_equal(length(blase_data@bins), 5)
    expect_equal(length(blase_data@pseudobulk_bins), 5)
})


test_that("generates pseudotime bins by pseudotime_range", {
    n_cells <- 150
    sce <- generate_test_sce(cells = n_cells, genes = 100)
    sce$pseudotime <- sqrt(seq_len(n_cells))

    blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 6, split_by = "pseudotime_range")

    cells_per_bin <- c(lapply(blase_data@pseudobulk_bins, ncol), recursive = TRUE)
    expect_equal(sum(cells_per_bin), n_cells)

    expected_output <- c(4, 12, 21, 29, 38, 46)
    expect_equal(cells_per_bin, expected_output)
})

test_that("generates pseudotime bins by cells", {
    n_cells <- 150
    sce <- generate_test_sce(cells = n_cells, genes = 100)
    sce$pseudotime <- sqrt(seq_len(n_cells))

    blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 6, split_by = "cells")

    cells_per_bin <- c(lapply(blase_data@pseudobulk_bins, ncol), recursive = TRUE)
    expect_equal(sum(cells_per_bin), n_cells)
    expected_output <- c(25, 25, 25, 25, 25, 25)
    expect_equal(cells_per_bin, expected_output)
})

test_that("throws error for invalid split_by parameter", {
    sce <- generate_test_sce(cells = 150, genes = 100)
    sce$pseudotime <- sqrt(seq_len(150))

    tmp1 <- function() {
        blase_data <- as.BlaseData(sce,
            pseudotime_slot = "pseudotime",
            n_bins = 10, split_by = "this_won't_work"
        )
    }
    expect_error(tmp1(), "split_by must be 'pseudotime_range' or 'cells'", fixed = TRUE)
})

test_that("throws error when pseudotime slot not available", {
    sce <- generate_test_sce(cells = 150, genes = 100)

    tmp1 <- function() {
        blase_data <- as.BlaseData(sce,
            pseudotime_slot = "pseudotime_not_here",
            n_bins = 10, split_by = "pseudotime_range"
        )
    }
    expect_error(tmp1(), "Pseudotime slot 'pseudotime_not_here' does not exist", fixed = TRUE)
})

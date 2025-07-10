test_that("get_top_n_genes with associationTest ignores genes with pvalue greater than cutoff", {
    assoRes <- data.frame(
        row.names = c("A", "B", "C", "D", "E"),
        waldStat = c(100, 50, 25, 10, 5),
        pvalue = c(0.01, 0.5, 0.005, 0.13, 0.039)
    )

    result <- get_top_n_genes(assoRes, n_genes = 3, p_cutoff = 0.04)

    expected_result <- c("A", "C", "E")
    expect_equal(length(result), length(expected_result))
    expect_equal(result, expected_result)
})

test_that("get_top_n_genes with associationTest should order genes by waldStat", {
    assoRes <- data.frame(
        row.names = c("A", "B", "C", "D"),
        waldStat = c(25, 50, 100, 10),
        pvalue = c(0.01, 0.5, 0.005, 0.13)
    )

    result <- get_top_n_genes(assoRes, n_genes = 2)

    expected_result <- c("C", "A")
    expect_equal(length(result), length(expected_result))
    expect_equal(result, expected_result)
})

test_that("get_top_n_genes with associationTest should provide less genes than
          asked for if fewer than that number meet the cutoff", {
    assoRes <- data.frame(
        row.names = c("A", "B", "C", "D"),
        waldStat = c(100, 50, 25, 10),
        pvalue = c(0.01, 0.5, 0.005, 0.13)
    )

    result <- get_top_n_genes(assoRes, n_genes = 10, p_cutoff = 0.8)
    expected_result <- c("A", "B", "C", "D")
    expect_equal(length(result), length(expected_result))
    expect_equal(result, expected_result)
})

test_that("get_bins_as_bulk will handle happy case", {
    sce <- generate_test_sce(cells = 500, genes = 50)
    sce$pseudotime <- seq_len(500)
    sce$pseudotime[1] <- 0
    sce <- assign_pseudotime_bins(sce, 5, pseudotime_slot = "pseudotime", split_by = "cells")
    sce$replicate <- rep(c(1, 2), 250)

    SingleCellExperiment::counts(sce[, sce$replicate == 2]) <- SingleCellExperiment::counts(sce[, sce$replicate == 2]) * 2

    result <- get_bins_as_bulk(sce, min_cells_for_bulk = 1, replicate_slot = "replicate")

    expected_result <- data.frame(
        bin_1_rep_1 = rep(150, 50),
        bin_1_rep_2 = rep(300, 50),
        bin_2_rep_1 = rep(150, 50),
        bin_2_rep_2 = rep(300, 50),
        bin_3_rep_1 = rep(150, 50),
        bin_3_rep_2 = rep(300, 50),
        bin_4_rep_1 = rep(150, 50),
        bin_4_rep_2 = rep(300, 50),
        bin_5_rep_1 = rep(150, 50),
        bin_5_rep_2 = rep(300, 50)
    )
    rownames(expected_result) <- paste0("G", rownames(expected_result))
    expected_result <- expected_result[order(rownames(expected_result)), ]

    expect_equal(result, expected_result)
})

test_that("get_bins_as_bulk will sample from replicates with small numbers of a given bin", {
    sce <- generate_test_sce(cells = 500, genes = 50)
    sce$pseudotime <- seq_len(500)
    sce$pseudotime[1] <- 0
    sce <- assign_pseudotime_bins(sce, 5, pseudotime_slot = "pseudotime", split_by = "cells")
    # The end of the trajectory is replicate 1 inflated
    sce$replicate <- c(rep(c(1, 2), 200), rep(1, 100))

    SingleCellExperiment::counts(sce[, sce$replicate == 2]) <- SingleCellExperiment::counts(sce[, sce$replicate == 2]) * 2

    result <- get_bins_as_bulk(sce, min_cells_for_bulk = 1, replicate_slot = "replicate")

    expected_result <- data.frame(
        bin_1_rep_1 = rep(150, 50),
        bin_1_rep_2 = rep(300, 50),
        bin_2_rep_1 = rep(150, 50),
        bin_2_rep_2 = rep(300, 50),
        bin_3_rep_1 = rep(150, 50),
        bin_3_rep_2 = rep(300, 50),
        bin_4_rep_1 = rep(150, 50),
        bin_4_rep_2 = rep(300, 50),
        bin_5_rep_1_1 = rep(225, 50),
        bin_5_rep_1_2 = rep(225, 50)
    )
    rownames(expected_result) <- paste0("G", rownames(expected_result))
    expected_result <- expected_result[order(rownames(expected_result)), ]

    expect_equal(result, expected_result)
})

test_that("get_bins_as_bulk will skip a bin if there aren't enough cells according to min_cells_for_bulk param", {
    result <- ""
    sce <- generate_test_sce(cells = 500, genes = 50)
    sce$pseudotime <- seq_len(500)
    sce$pseudotime[1] <- 0
    sce <- assign_pseudotime_bins(sce, 5, pseudotime_slot = "pseudotime", split_by = "cells")
    # The end of the trajectory is replicate 1 inflated
    sce$replicate <- c(rep(c(1, 2), 200), rep(c(1, 2, 3, 4, 5), 20))

    SingleCellExperiment::counts(sce[, sce$replicate == 2]) <- SingleCellExperiment::counts(sce[, sce$replicate == 2]) * 2

    tmp1 <- function() {
        return(get_bins_as_bulk(sce, min_cells_for_bulk = 30, replicate_slot = "replicate"))
    }
    expect_message(tmp1(), "Couldn't create pseudobulks due to too few cells in every replicate for bin 5")
    suppressMessages(result <- tmp1())

    expected_result <- data.frame(
        bin_1_rep_1 = rep(150, 50),
        bin_1_rep_2 = rep(300, 50),
        bin_2_rep_1 = rep(150, 50),
        bin_2_rep_2 = rep(300, 50),
        bin_3_rep_1 = rep(150, 50),
        bin_3_rep_2 = rep(300, 50),
        bin_4_rep_1 = rep(150, 50),
        bin_4_rep_2 = rep(300, 50)
    )

    # No bin 5 in this result
    rownames(expected_result) <- paste0("G", rownames(expected_result))
    expected_result <- expected_result[order(rownames(expected_result)), ]

    expect_equal(result, expected_result)
})

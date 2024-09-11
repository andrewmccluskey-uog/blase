## WIP Test select_genes_by_fourier_method
test_that("throws error when method not in 'power', 'amplitude', 'r2'", {
    sce <- generate_test_sce(cells = 40, genes = 5)
    sce$pseudotime <- (1:40) / 40

    SingleCellExperiment::counts(sce, withDimnames = FALSE) <- matrix(1:200, ncol = 40, nrow = 5) * 3
    SingleCellExperiment::normcounts(sce, withDimnames = FALSE) <- matrix(1:200, ncol = 40, nrow = 5)

    waves <- get_waves(
        sce = sce,
        pseudotime_slot = "pseudotime",
        n_cores = 1
    )

    tmp1 <- function() {
        select_genes_by_fourier_method(
            x = generate_test_sce(cells = 20, genes = 10),
            waves = waves,
            n_genes = 4,
            n_groups = 2,
            top_n_per_group = 1,
            method = "something_made_up",
            force_spread_selection = TRUE
        )
    }
    expect_error(tmp1(), "Requested method is not valid, must be one of ['power','amplitude','r2']", fixed = TRUE)
})

## Test get_waves
test_that("throws error when sce does not have pseudotime slot filled", {
    tmp1 <- function() {
        get_waves(
            sce = generate_test_sce(cells = 20, genes = 10),
            pseudotime_slot = "non_existent_slot",
            n_cores = 1
        )
    }
    expect_error(tmp1(), "Pseudotime slot 'non_existent_slot' does not exist", fixed = TRUE)
})

test_that("calculation of waves is correct", {
    sce <- generate_test_sce(cells = 40, genes = 5)
    sce$pseudotime <- (1:40) / 40

    SingleCellExperiment::counts(sce, withDimnames = FALSE) <- matrix(1:200, ncol = 40, nrow = 5) * 3
    SingleCellExperiment::normcounts(sce, withDimnames = FALSE) <- matrix(1:200, ncol = 40, nrow = 5)

    waves <- get_waves(
        sce = sce,
        pseudotime_slot = "pseudotime",
        n_cores = 1
    )

    expected_result <- data.frame(
        amplitude = c(63.7274742159119, 63.7274742159119, 63.7274742159119, 63.7274742159119, 63.7274742159119),
        phase = c(0.737499918678986, 0.737499918678986, 0.737499918678986, 0.737499918678986, 0.737499918678986),
        k = c(1, 1, 1, 1, 1),
        r2 = c(0.610705409013491, 0.610705409013491, 0.610705409013491, 0.610705409013491, 0.610705409013491),
        gene = c("1", "2", "3", "4", "5"),
        total_expression = c(3940, 3980, 4020, 4060, 4100),
        peak_expression = c(287, 289, 291, 293, 295),
        cellcount_in_peak = c(2, 2, 2, 2, 2),
        power = c(0.222046948487498, 0.220510291404539, 0.218994756755711, 0.217499911999699, 0.216025336325125)
    )
    rownames(expected_result) <- c("1", "2", "3", "4", "5")



    expect_equal(waves, expected_result)
})

## Test redim_matrix
test_that("a target size larger than original number of rows throws error", {
    tmp1 <- function() {
        redim_matrix(
            mat = matrix(1:9, nrow = 3, ncol = 3),
            target_height = 5,
            target_width = 2,
            n_core = 1
        )
    }
    expect_error(tmp1(), "Input matrix must be bigger than target width and height.", fixed = TRUE)
})

test_that("a target size larger than original number of columns throws error", {
    tmp1 <- function() {
        redim_matrix(
            mat = matrix(1:9, nrow = 3, ncol = 3),
            target_height = 2,
            target_width = 10,
            n_core = 1
        )
    }
    expect_error(tmp1(), "Input matrix must be bigger than target width and height.", fixed = TRUE)
})

test_that("scales matrix down", {
    new_matrix <- redim_matrix(
        mat = matrix(1:9, nrow = 3, ncol = 3),
        target_height = 2,
        target_width = 2,
        n_core = 1
    )

    expect_equal(new_matrix, matrix(c(3, 4, 6, 7), nrow = 2, ncol = 2))
})

## test map_all_best_bins
test_that("mapping all bulks runs without errors", {
    counts_matrix <- matrix(c(seq_len(120) / 10, seq_len(120) / 5), ncol = 48, nrow = 5)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(
        normcounts = counts_matrix, logcounts = log(counts_matrix)
    ))
    colnames(sce) <- seq_len(48)
    rownames(sce) <- as.character(seq_len(5))
    sce$cell_type <- c(rep("celltype_1", 24), rep("celltype_2", 24))

    sce$pseudotime <- seq_len(48)
    sce$pseudotime[1] <- 0
    blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 4)
    blase_data@genes <- as.character(seq_len(5))

    bulk_counts <- matrix(seq_len(15) * 10, ncol = 3, nrow = 5)
    colnames(bulk_counts) <- c("A", "B", "C")
    rownames(bulk_counts) <- as.character(seq_len(5))

    mappings <- map_all_best_bins(blase_data, bulk_counts, bootstrap_iterations = 1)

    expect_equal(length(mappings), 3)

    expect_equal(mappings[[1]]@bulk_name, "A")
    expect_equal(mappings[[2]]@bulk_name, "B")
    expect_equal(mappings[[3]]@bulk_name, "C")

    expect_equal(mappings[[1]]@best_bin, 1)
    expect_equal(mappings[[2]]@best_bin, 1)
    expect_equal(mappings[[3]]@best_bin, 1)
})



## WIP test map_best_bin

test_that("throws error if gene list is null", {
    sce <- generate_test_sce()
    blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 5)

    tmp1 <- function() map_best_bin(blase_data, "1", data.frame(), bootstrap_iterations = 1)
    expect_error(tmp1(), "No genes to map with. Please add something to the genes(blase_data) slot.", fixed = TRUE)
})

test_that("throws error if gene list has nothing matching the bulks", {
    sce <- generate_test_sce()
    blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 5)
    genes(blase_data) <- c("Something", "Made", "Up")

    tmp1 <- function() map_best_bin(blase_data, "1", data.frame(), bootstrap_iterations = 1)
    expect_error(tmp1(), "No genes in genes(blase_data) exist in the rows of the bulk dataframe, exiting.", fixed = TRUE)
})

test_that("mapping runs without errors", {
    counts_matrix <- matrix(c(seq_len(120) / 10, seq_len(120) / 5), ncol = 48, nrow = 5)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(
        normcounts = counts_matrix, logcounts = log(counts_matrix)
    ))
    colnames(sce) <- seq_len(48)
    rownames(sce) <- as.character(seq_len(5))
    sce$cell_type <- c(rep("celltype_1", 24), rep("celltype_2", 24))

    sce$pseudotime <- seq_len(48)
    sce$pseudotime[1] <- 0
    blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 4)
    blase_data@genes <- as.character(seq_len(5))

    bulk_counts <- matrix(seq_len(15) * 10, ncol = 3, nrow = 5)
    colnames(bulk_counts) <- c("A", "B", "C")
    rownames(bulk_counts) <- as.character(seq_len(5))

    mapping <- map_best_bin(blase_data, "A", bulk_counts, bootstrap_iterations = 1)

    expect_equal(mapping@best_bin, 1)
})

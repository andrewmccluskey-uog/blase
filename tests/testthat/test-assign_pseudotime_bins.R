# assign_pseudotime_bins (dataframe)
test_that("throws error as pseudotime bins cannot be applied to bulks", {
    bulk_counts <- matrix(seq_len(15) * 10, ncol = 3, nrow = 5)
    colnames(bulk_counts) <- c("A", "B", "C")
    rownames(bulk_counts) <- as.character(seq_len(5))
    bulk_counts <- as.data.frame(bulk_counts)

    tmp1 <- function() assign_pseudotime_bins(bulk_counts)
    expect_error(tmp1(), "Can't update bulk data, using each sample as bins.", fixed = TRUE)
})

# assign_pseudotime_bins (seurat)
test_that("adds pseudotime bins to a seurat object", {
    seurat <- generate_test_seurat()
    seurat <- assign_pseudotime_bins(seurat, "pseudotime_range", 2, pseudotime_slot = "pseudotime")
    expect_equal(unname(seurat$pseudotime_bin), c(rep(1, 50), rep(2, 50)))
})

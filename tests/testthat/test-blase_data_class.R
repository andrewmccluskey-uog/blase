## Tests show
test_that("blase_data_class show() prints what we expect", {
    cells <- 15
    genes <- 20

    sce <- generate_test_sce(cells, genes)
    blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 5)

    expect_output(
        print(blase_data),
        "Blase Data with:\n\tbins: c(1, 2, 3, 4, 5)\n\tselected genes: character(0)",
        fixed = TRUE
    )
})

test_that("blase_data_class show() prints genes if set", {
    cells <- 15
    genes <- 20

    sce <- generate_test_sce(cells, genes)
    blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 5)
    genes(blase_data) <- rownames(sce)[1:5]

    expect_output(
        print(blase_data),
        "Blase Data with:\n\tbins: c(1, 2, 3, 4, 5)\n\tselected genes: c(\"G1\", \"G2\", \"G3\", \"G4\", \"G5\")",
        fixed = TRUE
    )
})


## Test genes accessors
test_that("genes can be set using accessor", {
    cells <- 15
    genes <- 20

    sce <- generate_test_sce(cells, genes)
    blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 5)
    genes(blase_data) <- as.character(seq_len(genes))

    expect_equal(blase_data@genes, as.character(seq_len(genes)))
})

test_that("genes can be gotten using accessor", {
    cells <- 15
    genes <- 20

    sce <- generate_test_sce(cells, genes)
    blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 5)
    blase_data@genes <- as.character(seq_len(genes))

    expect_equal(genes(blase_data), as.character(seq_len(genes)))
})

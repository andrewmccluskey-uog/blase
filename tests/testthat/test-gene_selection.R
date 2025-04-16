generate_gene_selection_test_sce <- function() {
  ncells <- 70
  ngenes <- 20

  # Each gene should have mean around its gene number
  counts = c()
  for (i in seq_len(ngenes)) {
    counts = c(counts, dnorm(seq_len(ncells), mean=(ncells/i), sd=1))
  }

  counts_matrix <- matrix(
    counts,
    ncol = ncells,
    nrow = ngenes
  )
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(
    counts = counts_matrix * 3,
    normcounts = counts_matrix,
    logcounts = log(counts_matrix)
  ))
  colnames(sce) <- paste0("cell", seq_len(ncells))
  rownames(sce) <- paste0("gene", seq_len(ngenes))
  sce$cell_type <- c(
    rep("celltype_1", ncells / 2),
    rep("celltype_2", ncells / 2)
  )
  sce$pseudotime <- seq_len(ncells)
  return(sce)
}

# tests get_top_n_genes
test_that("gene_selection get_top_n_genes gets top n genes with no lineages", {
  assoRes <- data.frame(
    row.names = c("A", "B", "C", "D", "E"),
    waldStat = c(25, 50, 100, 10, 35),
    pvalue = c(0.01, 0.5, 0.005, 0.13, 0.02)
  )

  expect_equal(get_top_n_genes(assoRes, n_genes = 2, p_cutoff = 0.05),
               c("C", "E"))
})

test_that("gene_selection get_top_n_genes gets top n genes with a lineage", {
  assoRes <- data.frame(
    row.names = c("A", "B", "C", "D", "E"),
    waldStat_1 = c(25, 50, 100, 10, 35),
    pvalue_1 = c(0.01, 0.5, 0.005, 0.13, 0.02),
    waldStat_2 = c(12, 70, 35, 4, 100),
    pvalue_2 = c(0.9, 0.005, 0.5, 0.0013, 0.02)
  )

  expect_equal(get_top_n_genes(
    assoRes,
    n_genes = 2,
    p_cutoff = 0.05,
    lineage = 2
  ),
  c("E", "B"))
})

test_that("gene_selection get_top_n_genes gets top n genes with no lineages and adjusted p value cutoff", {
    assoRes <- data.frame(
      row.names = c("A", "B", "C", "D", "E"),
      waldStat = c(25, 50, 100, 10, 35),
      pvalue = c(0.01, 0.5, 0.005, 0.13, 0.02)
    )

    expect_equal(get_top_n_genes(assoRes, n_genes = 3, p_cutoff = 0.2),
                 c("C", "E", "A"))
  }
)

test_that("gene_selection get_top_n_genes gets top n genes when requested genes higher than actual results", {
    assoRes <- data.frame(
      row.names = c("A", "B", "C", "D", "E"),
      waldStat = c(25, 50, 100, 10, 35),
      pvalue = c(0.01, 0.5, 0.005, 0.13, 0.02)
    )

    expect_equal(get_top_n_genes(assoRes, n_genes = 10, p_cutoff = 0.1),
                 c("C", "E", "A"))
  }
)

# tests calculate_gene_peakedness
test_that("gene_selection calculate_gene_peakedness throws error if pseudotime slot not present", {
  sce = generate_gene_selection_test_sce()
  tmp1 <- function() calculate_gene_peakedness(sce, pseudotime_slot="nothing")
  expect_error(tmp1(), "Pseudotime slot not in object", fixed = TRUE)
})

test_that("gene_selection calculate_gene_peakedness calculates peakedness", {

  sce = generate_gene_selection_test_sce()

  gene_peakedness <- calculate_gene_peakedness(sce, pseudotime_slot="pseudotime")

  # we'll check the first ten
  gene_peakedness = gene_peakedness[1:10,]

  expect_equal(gene_peakedness$gene, paste0("gene", 1:10))
  expect_equal(gene_peakedness$peak_pseudotime, c(18.9, 18.9, 9.1, 70, 70, 0.7, 9.8, 11.9, 0.7, 0.7))
  expect_equal(gene_peakedness$mean_in_window, c(
    0.045635429, 0.053911890,
    0.061627168, 0,
    0, 0.000033458,
    0.050363657, 0.050296111,
    0.060492681, 0.099735570
  ), tolerance = 0.0001)
  expect_equal(gene_peakedness$mean_out_window, c(
    0.001154014, 0.005113787,
    0.012975915, 0.028203005,
    0.028324385, 0.021635033,
    0.013108241, 0.009857911,
    0.008691179, 0.006832208
  ), tolerance = 0.0001)
  expect_equal(gene_peakedness$ratio, c(
    39.544944523, 10.542458045,
    4.749350382, 0.0,
    0.0, 0.001546453,
    3.842136779, 5.102106562,
    6.960238789, 14.597852918
  ), tolerance = 0.0001)
  expect_equal(gene_peakedness$window_start, c(
    15.4, 15.4,
    5.6, 66.5,
    66.5, -2.8,
    6.3, 8.4,
    -2.8, -2.8
  ))
  expect_equal(gene_peakedness$window_end, c(
    22.4, 22.4,
    12.6, 73.5,
    73.5, 4.2,
    13.3, 15.4,
    4.2, 4.2
  ))
  expect_equal(gene_peakedness$deviance_explained, c(
    0.982064659, 0.833595976,
    0.718683962, 0.007329955,
    0.001096911, 0.001128665,
    0.016843770, 0.939128426,
    0.045942118, 0.044634839
  ), tolerance = 0.0001)

})

# tests smooth_gene
test_that("gene_selection smooth_gene throws error if pseudotime slot not present", {
  sce = generate_gene_selection_test_sce()
  gene_peakedness <- calculate_gene_peakedness(sce, pseudotime_slot="pseudotime")

  tmp1 <- function() smooth_gene(sce, "gene5", pseudotime_slot="nothing")
  expect_error(tmp1(), "Pseudotime slot not in object", fixed = TRUE)
})

test_that("gene_selection smooth_gene calculates gene smoothing results", {

  sce = generate_gene_selection_test_sce()
  gene_peakedness <- calculate_gene_peakedness(sce, pseudotime_slot="pseudotime")

  # smooth_gene
  smoothed_gene5 <- smooth_gene(
    sce, "gene5", pseudotime_slot = "pseudotime")

  smoothed_gene5 = as.vector(smoothed_gene5)
  expect_equal(smoothed_gene5[1], 0.022075948, ignore_attr=TRUE)
  expect_equal(smoothed_gene5[81], 0.029716046, ignore_attr=TRUE)

})

# tests plot_gene_peakedness
test_that("gene_selection plot_gene_peakedness throws error if pseudotime slot not present", {
  sce = generate_gene_selection_test_sce()
  gene_peakedness <- calculate_gene_peakedness(sce, pseudotime_slot="pseudotime")

  tmp1 <- function() plot_gene_peakedness(sce, gene_peakedness, "gene5", pseudotime_slot="nothing")
  expect_error(tmp1(), "Pseudotime slot not in object", fixed = TRUE)
})

test_that("gene_selection plot_gene_peakedness throws error if gene not found", {
  sce = generate_gene_selection_test_sce()
  gene_peakedness <- calculate_gene_peakedness(sce, pseudotime_slot="pseudotime")

  tmp1 <- function() plot_gene_peakedness(sce, gene_peakedness, "made_up_gene", pseudotime_slot="pseudotime")
  expect_error(tmp1(), "Gene not in gene_peakedness_df, please make sure gene exists in dataset.", fixed = TRUE)
})

test_that("gene_selection plot_gene_peakedness throws error if gene in gene_peakedness_df more than once", {
  sce = generate_gene_selection_test_sce()
  gene_peakedness <- calculate_gene_peakedness(sce, pseudotime_slot="pseudotime")
  gene_peakedness <- rbind(gene_peakedness, gene_peakedness[1,])

  tmp1 <- function() plot_gene_peakedness(sce, gene_peakedness, gene_peakedness$gene[1], pseudotime_slot="pseudotime")
  expect_error(tmp1(), "Multiple copies of gene in gene_peakedness_df, please make sure only one exists.", fixed = TRUE)
})

test_that("gene_selection plot_gene_peakedness runs without error", {

  sce = generate_gene_selection_test_sce()
  gene_peakedness <- calculate_gene_peakedness(sce, pseudotime_slot="pseudotime")

  tmp1 = function() plot_gene_peakedness(sce, gene_peakedness, "gene5", pseudotime_slot="pseudotime")
  expect_no_error(tmp1)

})

# gene_peakedness_spread_selection
test_that("gene_selection gene_peakedness_spread_selection selects genes", {
  sce = generate_gene_selection_test_sce()
  gene_peakedness <- calculate_gene_peakedness(sce, pseudotime_slot="pseudotime")

  genes_to_use <- gene_peakedness_spread_selection(
    sce, gene_peakedness, genes_per_bin = 2, n_gene_bins = 1,
    pseudotime_slot="pseudotime")

  expect_equal(genes_to_use, c("gene11", "gene19"))

})

test_that("gene_selection gene_peakedness_spread_selection selects genes", {
  sce = generate_gene_selection_test_sce()
  gene_peakedness <- calculate_gene_peakedness(sce, pseudotime_slot="pseudotime")

  genes_to_use <- gene_peakedness_spread_selection(
    sce, gene_peakedness, genes_per_bin = 1, n_gene_bins = 2,
    pseudotime_slot="pseudotime")

  expect_equal(genes_to_use, c("gene11", "gene16"))
})

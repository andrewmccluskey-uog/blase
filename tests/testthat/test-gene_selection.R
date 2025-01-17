generate_gene_selection_test_sce <- function() {
  ncells <- 70
  ngenes <- 10
  counts_matrix <- matrix(
    c(seq_len(350) / 10, seq_len(350) / 5),
    ncol = ncells,
    nrow = ngenes
  )
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(
    normcounts = counts_matrix, logcounts = log(counts_matrix)
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

  expect_equal(gene_peakedness$gene, paste0("gene", 1:10))
  expect_equal(gene_peakedness$peak_pseudotime, rep_len(70, 10))
  expect_equal(gene_peakedness$mean_in_window, c(65.2, 65.4, 65.6, 65.8, 66, 66.2, 66.4, 66.6, 66.8, 67))
  expect_equal(gene_peakedness$mean_out_window, c(
    23.25303030, 23.40000000,
    23.54696970, 23.69393939,
    23.84090909, 23.98787879,
    24.13484848, 24.28181818,
    24.42878788, 24.57575758
  ))
  expect_equal(gene_peakedness$ratio, c(
    2.80393562, 2.79487179,
    2.78592111, 2.77708147,
    2.76835081, 2.75972713,
    2.75120849, 2.74279296,
    2.73447870, 2.72626387
  ))
  expect_equal(gene_peakedness$window_start, rep_len(66.5, 10))
  expect_equal(gene_peakedness$window_end, rep_len(73.5, 10))
  expect_equal(gene_peakedness$deviance_explained, c(
    0.92270963, 0.92458661,
    0.92647701, 0.92838311,
    0.93030661, 0.93180178,
    0.93324144, 0.93463523,
    0.93599024, 0.93731242
  ))

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

  expected_result = c(
    1.4912445963909,
    1.75244505491791,
    2.05837028285099,
    2.41446414496174,
    2.82595943037959,
    3.29752981626186,
    3.83282443273115,
    4.4339036265525,
    5.10059030161558,
    5.82976759451327)

  expect_equal(as.vector(smoothed_gene5[1:10]), expected_result)

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

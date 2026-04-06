# Evaluate n_bins and n_genes for bin mapping

Will use the n_bins and n_genes implied by the `sce` and
`pseudotime_bins_top_n_genes_df` parameters and return quality metrics
and an optional chart.

## Usage

``` r
evaluate_parameters(
  blase_data,
  bootstrap_iterations = 200,
  BPPARAM = BiocParallel::SerialParam(),
  make_plot = FALSE,
  plot_columns = 4
)
```

## Arguments

- blase_data:

  The
  [BlaseData](https://andrewmccluskey-uog.github.io/blase/reference/BlaseData-class.md)
  object to use.

- bootstrap_iterations:

  Integer. Iterations for bootstrapping when calculating strong
  mappings.

- BPPARAM:

  The
  [BiocParallel::BiocParallelParam](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  configuration. Defaults to
  [BiocParallel::SerialParam](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)

- make_plot:

  Boolean. Whether or not to render the plot showing the correlations
  for each pseudobulk bin when we try to map the given bin.

- plot_columns:

  Integer. How many columns to use in the plot.

## Value

A vector of length 3:

- "worst top 2 distance" decimal containing the lowest difference
  between the absolute values of the top 2 most correlated bins for each
  bin. Higher is better for differentiating.

- "mean top 2 distance" decimal containing the mean top 2 distance
  across the entire set of genes and bins. Higher is better for
  differentiation, but it should matter less than the worst value.

- "strong_mapping_pct" decimal from 0-1. The percent of mappings for
  this setup which were annotated as strong by BLASE.

## Examples

``` r
ncells <- 70
ngenes <- 100
counts_matrix <- matrix(
    c(seq_len(3500) / 10, seq_len(3500) / 5),
    ncol = ncells,
    nrow = ngenes
)
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(
    normcounts = counts_matrix, logcounts = log(counts_matrix)
))
colnames(sce) <- seq_len(ncells)
rownames(sce) <- as.character(seq_len(ngenes))
sce$cell_type <- c(
    rep("celltype_1", ncells / 2),
    rep("celltype_2", ncells / 2)
)

sce$pseudotime <- seq_len(ncells) - 1
genelist <- as.character(seq_len(ngenes))

# Evaluating created BlaseData
blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 10)
genes(blase_data) <- genelist[1:20]

# Check convexity of parameters
evaluate_parameters(blase_data, make_plot = TRUE)
#> Warning: Bulk ID matches a gene, if this fails then check you areusing bulk name and not geneIds:1
#> Warning: Genes for mapping not all in bulk, using 0 genes available in both reference and bulk.
#> Error: BiocParallel errors
#>   1 remote errors, element index: 1
#>   9 unevaluated and other errors
#>   first remote error:
#> Error in if (ncol(pseudobulk_bins(blase_data)[[i]]) <= 1) {: argument is of length zero
```

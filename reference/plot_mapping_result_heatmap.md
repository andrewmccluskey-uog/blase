# Plot a mapping result heatmap

Plots Spearman's Rho as the fill colour, and adds \* if the
[MappingResult](https://andrewmccluskey-uog.github.io/blase/reference/MappingResult.md)
was strongly assigned.

## Usage

``` r
plot_mapping_result_heatmap(
  mapping_result_list,
  heatmap_fill_scale = NULL,
  annotate_strong = TRUE,
  annotate_correlation = FALSE,
  bin_order = NULL,
  text_background = FALSE
)
```

## Arguments

- mapping_result_list:

  A list of
  [MappingResult](https://andrewmccluskey-uog.github.io/blase/reference/MappingResult.md)
  objects to include in the heatmap.

- heatmap_fill_scale:

  The ggplot2 compatible fill gradient scale to apply to the heatmap.

- annotate_strong:

  Boolan. Whether to annotate the heatmap with strong results or not,
  defaults to TRUE.

- annotate_correlation:

  Boolean. Whether to annotate the heatmap with the correlation of bin
  to each bulk sample. Defaults to FALSE.

- bin_order:

  Vector of integers. A vector of the bin ids in which to plot the
  pseudotime bins along the x-axis.

- text_background:

  Boolean. Whether to show background on labels or not. Has no effect if
  no annotations are enabled.

## Value

A
[ggplot2::ggplot2](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
heatmap showing the correlations of each mapping result across every
pseudotime bin.

## Examples

``` r
counts_matrix <- matrix(
    c(seq_len(120) / 10, seq_len(120) / 5),
    ncol = 48, nrow = 5
)
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(
    normcounts = counts_matrix, logcounts = log(counts_matrix)
))
colnames(sce) <- seq_len(48)
rownames(sce) <- as.character(seq_len(5))
sce$cell_type <- c(rep("celltype_1", 24), rep("celltype_2", 24))

sce$pseudotime <- seq_len(48) - 1
blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 4)
genes(blase_data) <- as.character(seq_len(5))

bulk_counts <- matrix(seq_len(15) * 10, ncol = 3, nrow = 5)
colnames(bulk_counts) <- c("A", "B", "C")
rownames(bulk_counts) <- as.character(seq_len(5))

# Map to bin
result <- map_best_bin(blase_data, "B", bulk_counts)
result
#> MappingResult for 'B': best_bin=1 correlation=1 top_2_distance=0
#>   Strong Result: FALSE (next max upper  1 )
#>   with history for scores against 4  bins
#>   Bootstrapped with 200 iterations

# Map all bulks to bin
results <- map_all_best_bins(blase_data, bulk_counts)

# Plot Heatmap
plot_mapping_result_heatmap(list(result))
#> Inferred correlation metric.


# Plot Correlation
plot_mapping_result_corr(result)


# Plot populations
sce <- assign_pseudotime_bins(
    sce,
    pseudotime_slot = "pseudotime", n_bins = 4
)
plot_bin_population(sce, best_bin(result), group_by_slot = "cell_type")


# Getters
bulk_name(result)
#> [1] "B"
best_bin(result)
#> [1] 1
best_correlation(result)
#> [1] 1
top_2_distance(result)
#> [1] 0
strong_mapping(result)
#> [1] FALSE
mapping_history(result)
#>   bin correlation lower_bound upper_bound
#> 1   1           1           1           1
#> 2   2           1           1           1
#> 3   3           1           1           1
#> 4   4           1           1           1
bootstrap_iterations(result)
#> [1] 200
metric(result)
#> [1] "spearman"

# Setters
bulk_name(result) <- "New Name"
```

# Annotate a SCE with BLASE Mappings

Annotates an SCE with the names of bulk samples that best match each
pseudotime bin. For each pseudotime bin, we find the highest correlation
with a bulk sample that was mapped against it. Because of this approach,
a bulk which mapped best to another pseudotime bin may be the best
correlation with the current pseudotime bin of interest.

## Usage

``` r
annotate_sce(
  sce,
  blase_results,
  annotation_col = "BLASE_Annotation",
  include_stats = FALSE
)
```

## Arguments

- sce:

  The
  [SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  to annotate.

- blase_results:

  A list of
  [MappingResult](https://andrewmccluskey-uog.github.io/blase/reference/MappingResult.md)
  to use for the annotation.

- annotation_col:

  String. The name of the metadata column in which to store the new
  annotations.

- include_stats:

  Boolean. Whether or not to include metadata columns containing the
  correlation of the best matching bin, and whether that mapping was
  strong.

## Value

A
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
with annotations added to metadata (in a column defined by
`annotation_col`), and the correlations in
`BLASE_Annotation_Correlation` if `include_stats` is enabled.

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

# Map all bulks to bin
results <- map_all_best_bins(blase_data, bulk_counts)

sce <- assign_pseudotime_bins(
    sce,
    pseudotime_slot = "pseudotime", n_bins = 4
)

# Annotate SC from existing bulk
sce <- annotate_sce(sce, results)
table(sce$BLASE_Annotation)
#> 
#>  A 
#> 48 
```

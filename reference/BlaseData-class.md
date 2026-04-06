# Blase Data Object

For creation details, see
[`as.BlaseData()`](https://andrewmccluskey-uog.github.io/blase/reference/as.BlaseData.md)

## Value

A BlaseData object

## Slots

- `pseudobulk_bins`:

  list of [data.frame](https://rdrr.io/r/base/data.frame.html)s. Each
  item is a normalised count matrix representing a bin, where a column
  is a cell in the bin and each row is a gene.

- `bins`:

  list. A list of bin names for each timepoint.

- `genes`:

  list. A list of the genes selected for discriminating timepoints.

## Examples

``` r
counts <- matrix(rpois(100, lambda = 10), ncol = 10, nrow = 10)
sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(normcounts = counts)
)
sce$pseudotime <- seq_len(10) - 1
data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 3)
genes(data) <- as.character(seq_len(10))

genes(data)
#>  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10"
```

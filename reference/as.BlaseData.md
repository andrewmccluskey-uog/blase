# Conversion to BlaseData

Conversion to BlaseData

## Usage

``` r
as.BlaseData(x, ...)

# S4 method for class 'SingleCellExperiment'
as.BlaseData(
  x,
  pseudotime_slot = "slingPseudotime_1",
  n_bins = 20,
  split_by = "pseudotime_range"
)
```

## Arguments

- x:

  An object to take counts from

- ...:

  additional arguments passed to object-specific methods.

- pseudotime_slot:

  String or vector of strings. The
  [SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  slot(s) containing pseudotime values for each cell to be passed to
  [`assign_pseudotime_bins()`](https://andrewmccluskey-uog.github.io/blase/reference/assign_pseudotime_bins.md).

- n_bins:

  Integer. The number of bins to create, passed to
  [`assign_pseudotime_bins()`](https://andrewmccluskey-uog.github.io/blase/reference/assign_pseudotime_bins.md).

- split_by:

  String. The split_by method to be passed on to
  [`assign_pseudotime_bins()`](https://andrewmccluskey-uog.github.io/blase/reference/assign_pseudotime_bins.md).
  Must be one of `pseudotime_range` or `cells`.

## Value

An
[BlaseData](https://andrewmccluskey-uog.github.io/blase/reference/BlaseData-class.md)
object

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

# Show an BlaseData object

Show an BlaseData object

## Usage

``` r
# S4 method for class 'BlaseData'
show(object)
```

## Arguments

- object:

  a
  [BlaseData](https://andrewmccluskey-uog.github.io/blase/reference/BlaseData-class.md)
  object

## Value

A character vector describing the BLASE object

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

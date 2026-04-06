# Get pseudobulk bins of a BLASE Data object.

Get pseudobulk bins of a BLASE Data object.

## Usage

``` r
pseudobulk_bins(x)

# S4 method for class 'BlaseData'
pseudobulk_bins(x)
```

## Arguments

- x:

  a
  [BlaseData](https://andrewmccluskey-uog.github.io/blase/reference/BlaseData-class.md)
  object

## Value

List of dataframes. Each dataframe is the normalised counts of cells by
genes in each pseudotime bin. List index is the pseudotime bin.

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

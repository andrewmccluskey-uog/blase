# Set name of bulk of a BLASE Mapping Results object.

Set name of bulk of a BLASE Mapping Results object.

## Usage

``` r
bulk_name(x) <- value

# S4 method for class 'MappingResult'
bulk_name(x) <- value
```

## Arguments

- x:

  a
  [MappingResult](https://andrewmccluskey-uog.github.io/blase/reference/MappingResult.md)
  object

- value:

  String. The name of the bulk used to map against.

## Value

Nothing

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

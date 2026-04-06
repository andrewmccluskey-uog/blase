# Set genes of a BLASE Data object.

Set genes of a BLASE Data object.

## Usage

``` r
pseudobulk_bins(x) <- value

# S4 method for class 'BlaseData'
pseudobulk_bins(x) <- value
```

## Arguments

- x:

  a
  [BlaseData](https://andrewmccluskey-uog.github.io/blase/reference/BlaseData-class.md)
  object

- value:

  List of dataframes. Each dataframe is the normalised counts of cells
  by genes in each pseudotime bin. List index is the pseudotime bin.

## Value

Nothing

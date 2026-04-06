# Get a pseudobulk of bins with at least 2 replicates

This function will try to create a pseudobulked count matrix for the
bins. When a replicate has too few cells, it is discounted. If only one
exists, then we sample from it twice to create the pseudobulks.

## Usage

``` r
get_bins_as_bulk(
  pseudotime_sce,
  min_cells_for_bulk = 50,
  replicate_slot = "replicate"
)
```

## Arguments

- pseudotime_sce:

  The
  [SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object to get the bins from

- min_cells_for_bulk:

  Integer. The minimum cells to look for per replicate and bin.

- replicate_slot:

  String. The name of the matadata column in the Single Cell Experiment
  that contains replicate information

## Value

A dataframe containing the pseudobulked counts matrix.

## Examples

``` r
library(SingleCellExperiment, quietly = TRUE)
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
library(blase)
counts <- matrix(rpois(1000, lambda = 10), ncol = 100, nrow = 10)
sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(normcounts = counts, counts = counts / 2)
)
sce$pseudotime <- seq_len(100) - 1
colnames(sce) <- seq_len(100)
rownames(sce) <- as.character(seq_len(10))
sce <- assign_pseudotime_bins(sce,
    n_bins = 5,
    pseudotime_slot = "pseudotime", split_by = "cells"
)
sce$replicate <- rep(c(1, 2), 50)
result <- get_bins_as_bulk(
    sce,
    min_cells_for_bulk = 1,
    replicate_slot = "replicate"
)
result
#>    bin_1_rep_1 bin_1_rep_2 bin_2_rep_1 bin_2_rep_2 bin_3_rep_1 bin_3_rep_2
#> 1         49.5        53.0        58.0        58.5        48.5        53.0
#> 10        53.0        47.0        49.5        44.5        48.5        45.0
#> 2         53.0        58.5        50.0        60.0        52.5        52.5
#> 3         51.0        49.0        56.5        41.5        49.0        53.0
#> 4         55.0        49.5        49.0        47.0        55.0        45.5
#> 5         52.5        45.0        49.0        51.0        51.5        54.5
#> 6         42.0        45.0        54.5        54.0        48.5        36.5
#> 7         50.0        49.0        44.5        45.0        46.0        48.5
#> 8         50.0        52.0        43.5        54.0        56.5        52.5
#> 9         48.0        48.0        50.5        56.5        50.0        57.5
#>    bin_4_rep_1 bin_4_rep_2 bin_5_rep_1 bin_5_rep_2
#> 1         50.5        51.0        48.5        46.0
#> 10        53.5        53.0        55.0        44.5
#> 2         47.0        49.5        43.0        52.0
#> 3         52.0        49.5        40.0        47.0
#> 4         51.0        51.0        46.5        45.0
#> 5         57.0        54.0        48.0        54.0
#> 6         53.0        45.5        44.5        52.0
#> 7         47.0        54.5        44.0        47.5
#> 8         49.5        38.0        57.0        44.0
#> 9         44.5        51.5        54.0        52.5
```

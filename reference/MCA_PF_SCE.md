# Malaria Cell Atlas Plasmodium falciparum for BLASE Vignette

Data from the Malaria Cell Atlas, with the following additional
processing:

1.  Genes renamed to match bulk samples in vignette

2.  Subset to 2500 cells

3.  Normalised

4.  Highly variable genes identified

5.  Pseudotime calculated

6.  Genes subset to include a spread of those found to have high ratios
    by BLASE's "Gene Peakedness" measure.

## Usage

``` r
MCA_PF_SCE
```

## Format

An object of class `SingleCellExperiment` with 1746 rows and 2500
columns.

## Source

<https://www.malariacellatlas.org/atlas/plasmodium-falciparum-atlas/>

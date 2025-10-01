#' TradeSeq Example SCE for BLASE Vignette
#'
#' Data from the TradeSeq vignette, with the following additional
#' processing applied:
#'
#' 1. Pseudotime calculated
#' 2. TradeSeq applied
#' 3. Log normalised and normalised counts calculated
#' 4. Erythrocyte cell type removed
#' 5. UMAP calculated
#'
#' @concept data
#' @source <https://bioconductor.org/packages/devel/bioc/vignettes/tradeSeq/inst/doc/tradeSeq.html>
"tradeSeq_BLASE_example_sce"

#' Malaria Cell Atlas Plasmodium falciparum for BLASE Vignette
#'
#' Data from the Malaria Cell Atlas, with the following additional processing:
#' 1. Genes renamed to match bulk samples in vignette
#' 2. Subset to 2500 cells
#' 3. Normalised
#' 4. Highly variable genes identified
#' 5. Pseudotime calculated
#' 6. Genes subset to include a spread of those found to have
#'    high ratios by BLASE's "Gene Peakedness" measure.
#'
#' @concept data
#' @source <https://www.malariacellatlas.org/atlas/plasmodium-falciparum-atlas/>
"MCA_PF_SCE"

#' Painter 2018 Plasmodium falciparum 48h asexual lifecycle microarray data
#'
#' Data originally from <https://doi.org/10.1038/s41467-018-04966-3>. Used as
#' generated in the BLASE reproducibility documents available at
#' <https://zenodo.org/records/16615703>, however genes have been subset to
#' reduce file size.
#'
#' @concept data
#' @source <https://zenodo.org/records/16615703>
"painter_microarray"

#' Zhang 2021 Plasmodium falciparum heat shock bulk data
#'
#' Data originally from <https://doi.org/10.1038/s41467-021-24814-1>. Used as
#' generated in the BLASE reproducibility documents available at
#' <https://zenodo.org/records/16615703>, however genes have been subset to
#' reduce file size.
#'
#' @concept data
#' @source <https://zenodo.org/records/16615703>
"zhang_2021_heat_shock_bulk"

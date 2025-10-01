## code to prepare `painter_microarray` dataset goes here

# Download and load the "blase_reproducibility/R/data/painter_2018_pf.rds
# from https://doi.org/10.5281/zenodo.16615703

painter_microarray <- readRDS("vignettes/articles/Data/painter_2018_pf.rds")
data(MCA_PF_SCE)
painter_microarray = painter_microarray[rownames(painter_microarray) %in% rownames(MCA_PF_SCE),]

usethis::use_data(painter_microarray, overwrite = TRUE, compress = "xz")

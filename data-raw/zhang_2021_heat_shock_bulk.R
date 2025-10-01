## code to prepare `zhang_2021_heat_shock_bulk` dataset goes here

# Download and load the
# "blase_reproducibility/R/data/zhang_2021_heat_shock_bulk.rds"
# from https://doi.org/10.5281/zenodo.16615703

zhang_2021_heat_shock_bulk <- readRDS("zhang_2021_heat_shock_bulk.rds")
data(MCA_PF_SCE)
zhang_2021_heat_shock_bulk = zhang_2021_heat_shock_bulk[rownames(zhang_2021_heat_shock_bulk) %in% rownames(MCA_PF_SCE),]

usethis::use_data(zhang_2021_heat_shock_bulk, overwrite = TRUE, compress = "xz")

## code to prepare `MCA_PF_SCE` dataset goes here

library(scran)
library(slingshot)
library(BiocParallel)
library(fs)
library(patchwork)
library(scater)
library(utils)
library(ami)
library(blase)

RNGversion("3.5.0")
SEED <- 7
set.seed(SEED)
N_CORES <- 4

mca_path <- fs::path("MCA", ext = "zip")
if (!file.exists(mca_path)) {
  download.file("https://www.malariacellatlas.org/downloads/pf.zip", mca_path)
  unzip(mca_path, exdir = paste0("MCA"))
} else {
  print("Using cached")
}

mca_counts_path <- fs::path("MCA", "pf-ch10x-raw", ext = "csv")
sc_readcounts <- read.csv(mca_counts_path, row.names = 1)

mca_annotation_path <- fs::path("MCA", "pf-ch10x-data", ext = "csv")
sc_annotations <- read.csv(mca_annotation_path, row.names = 1)

rownames(sc_annotations) <- gsub(x = rownames(sc_annotations), pattern = "-", replacement = ".", fixed = TRUE)

genes_to_fix <- rownames(sc_readcounts)[!(rownames(sc_readcounts) %in% sub(x = rownames(sc_readcounts), pattern = "\\.[0-9]", replacement = ""))]
genes_to_fix <- unique(sub(x = genes_to_fix, pattern = "\\.[0-9]", replacement = ""))

new_rows <- data.frame()
rownames_to_remove <- c()

# Do some pre-work to get colSums quickly
sc_readcounts_matrix <- data.matrix(sc_readcounts)
sc_readcounts_n <- ncol(sc_readcounts_matrix)

for (gene in genes_to_fix) {
  targetRowNames <- rownames(sc_readcounts)[grep(x = rownames(sc_readcounts), pattern = gene)]
  rownames_to_remove <- c(rownames_to_remove, targetRowNames)
  # N = col index, M = row index
  counts <- .colSums(sc_readcounts_matrix[targetRowNames, ], m = length(targetRowNames), n = sc_readcounts_n)
  counts <- t(data.frame(counts = counts))
  rownames(counts) <- c(gene)
  colnames(counts) <- colnames(sc_readcounts)
  new_rows <- rbind(new_rows, counts)
}

sc_readcounts <- sc_readcounts[!(rownames(sc_readcounts) %in% rownames_to_remove), ]
sc_readcounts <- rbind(sc_readcounts, new_rows)

genes_to_fix <- rownames(sc_readcounts)[!(rownames(sc_readcounts) %in% sub(x = rownames(sc_readcounts), pattern = "\\.[0-9]", replacement = ""))]
genes_to_fix <- unique(sub(x = genes_to_fix, pattern = "\\.[0-9]", replacement = ""))
print(paste("Remaining genes with dups:", genes_to_fix))

rm(new_rows, rownames_to_remove, sc_readcounts_matrix, sc_readcounts_n, genes_to_fix, targetRowNames, gene)
gc()

sparse_matrix <- as(as.matrix(sc_readcounts), "sparseMatrix")
sce <- SingleCellExperiment(assays = list(counts = sparse_matrix), colData = sc_annotations)
rownames(sce) <- sub(x = rownames(sce), pattern = "\\.[0-9]", replacement = "")

colData(sce)$STAGE_LR <- factor(colData(sce)$STAGE_LR)
colData(sce)$STAGE_HR <- factor(colData(sce)$STAGE_HR)
colData(sce)$DAY <- factor(colData(sce)$DAY)
colData(sce)$STRAIN <- NULL
colData(sce)$STAGE_HR2 <- factor(colData(sce)$STAGE_HR2)
colData(sce)$HOST <- NULL
colData(sce)$CLUSTER <- factor(colData(sce)$CLUSTER)

sce <- subset(sce, , STAGE_LR != "gametocyte")

# Subsample to reduce memory consumption
cells_to_keep <- sample(colnames(sce), size = 2500)
sce <- sce[, colnames(sce) %in% cells_to_keep]
rm(sc_annotations, sc_readcounts, cells_to_keep)

sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)
normcounts(sce) <- exp(logcounts(sce))

dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop = 0.3)

reducedDim(sce, "PCA") <- as.matrix(cbind(colData(sce)["PC_1"], colData(sce)["PC_2"], colData(sce)["PC_3"]))
# better separation in the 2nd and 3rd dimension for just asexual stages
reducedDim(sce, "UMAP") <- as.matrix(cbind(colData(sce)["UMAP_2"], colData(sce)["UMAP_3"]))

sce <- slingshot(sce, reducedDim = "UMAP", clusterLabels = "STAGE_HR", start.clus = "early ring")

gene_peakedness_info <- calculate_gene_peakedness(
  sce,
  window_pct = 5,
  knots = 18,
  BPPARAM = SerialParam()
)

genes_to_keep <- gene_peakedness_spread_selection(
  sce,
  gene_peakedness_info,
  genes_per_bin = 65,
  n_gene_bins = 30
)

# We remove some genes here to reduce file size.
sce <- sce[genes_to_keep, ]

logcounts(sce) <- NULL
sce$slingshot <- NULL

MCA_PF_SCE<-sce

usethis::use_data(MCA_PF_SCE, overwrite = TRUE, compress = "xz")

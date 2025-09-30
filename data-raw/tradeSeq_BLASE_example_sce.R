## code to prepare `tradeSeq_BLASE_example_sce` dataset goes here
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")

pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
tradeSeq_BLASE_example_sce <- fitGAM(
  counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
  nknots = 6, verbose = FALSE, BPPARAM = MulticoreParam(N_CORES)
)
tradeSeq_BLASE_example_sce$pseudotime <- pseudotime[, "curve1"]
tradeSeq_BLASE_example_sce$celltype <- celltype

tradeSeq_BLASE_example_sce <- subset(tradeSeq_BLASE_example_sce, , celltype != "Erythrocyte")

tradeSeq_BLASE_example_sce <- computeSumFactors(tradeSeq_BLASE_example_sce)
tradeSeq_BLASE_example_sce <- logNormCounts(tradeSeq_BLASE_example_sce)
normcounts(tradeSeq_BLASE_example_sce) <- exp(logcounts(tradeSeq_BLASE_example_sce))
tradeSeq_BLASE_example_sce <- runUMAP(tradeSeq_BLASE_example_sce)

usethis::use_data(tradeSeq_BLASE_example_sce, overwrite = TRUE)

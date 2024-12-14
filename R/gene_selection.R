#' Get Top Genes From An AssociationTestResult
#'
#' Pulls the genes with the highest wald statistic from an association
#' test result, with a p value cutoff.
#'
#' @concept gene-selection
#'
#' @param association_test_results The association test results data frame
#' to take the genes from.
#' @param n_genes The number of genes to return. Defaults to 40.
#' @param lineage The Lineage to use. The Defaults to NA, which assumes
#' the test was run with Lineages=False.
#' @param p_cutoff The P value cutoff to use. Defaults to less than 0.05.
#'
#' @return A vector of the names of the genes that best describe
#' a lineage's trajectory.
#' @export
#'
#' @examples
#' assoRes <- data.frame(
#'     row.names = c("A", "B", "C", "D"),
#'     waldStat = c(25, 50, 100, 10),
#'     pvalue = c(0.01, 0.5, 0.005, 0.13)
#' )
#' get_top_n_genes(assoRes, n_genes = 2)
get_top_n_genes <- function(
    association_test_results,
    n_genes = 40,
    lineage = NA,
    p_cutoff = 0.05) {
  pvalue_slot_for_lineage <- "pvalue"
  wald_slot_for_lineage <- "waldStat"

  if (!is.na(lineage)) {
    pvalue_slot_for_lineage <- paste0("pvalue_", lineage)
    wald_slot_for_lineage <- paste0("waldStat_", lineage)
  }

  # P Cutoff
  asso_results_copy <- association_test_results[
    association_test_results[, pvalue_slot_for_lineage] < p_cutoff,
  ]
  # Remove NAs
  asso_results_copy <- asso_results_copy[
    !is.na(rownames(asso_results_copy)),
  ]
  # Sort by wald stat
  asso_results_copy <- asso_results_copy[
    order(-asso_results_copy[, wald_slot_for_lineage]),
  ]

  topPdtGenesNames <- rownames(asso_results_copy)[seq_len(n_genes)]
  topPdtGenesNames <- topPdtGenesNames[!is.na(topPdtGenesNames)]
  return(topPdtGenesNames)
}


#' calculate_gene_peakedness
#'
#' @description
#' Calculate the peakedness of a gene. The power is the ratio of the mean of reads
#' 5% either side of the smoothed peak of the gene's expression over pseudotime
#' against the mean of the reads outside of this.
#'
#' @param sce SCE to do the calculations on.
#' @param window_pct the size of the window to consider, as a percentage
#' of the maximum pseudotime value.
#' @param BPPARAM The BiocParallel parameter for parallelisation
#' Defaults to [BiocParallel::SerialParam].
#'
#' @return Dataframe, where each row is a gene, and the following columns:
#' mean_expression_in_window, mean_expression_out_window, ratio
#' @export
#'
#' @concept gene-selection
#' @import BiocParallel
#' @import mgcv
#'
#' @examples
calculate_gene_peakedness <- function(
    sce,
    window_pct = 10,
    pseudotime_slot="slingPseudotime_1",
    knots = 10,
    BPPARAM = BiocParallel::SerialParam()) {

  if (!(pseudotime_slot %in% colnames(colData(sce)))) {
    stop("Pseudotime slot not in object")
  }

  pseudotime <- colData(sce)[[pseudotime_slot]]
  normalised_counts = normcounts(sce)

  dataframes = list()
  for (gene_id in rownames(normalised_counts)) {
    df = as.data.frame(normalised_counts[gene_id,])
    colnames(df) <- c(gene_id)
    rownames(df) <- colnames(normalised_counts)
    dataframes[[length(dataframes)+1]] = df
  }

  results <- BiocParallel::bplapply(
    dataframes,
    function(df) {

      gene = colnames(df)[1]
      to_smooth = data.frame(nc=df[,gene], pdt=pseudotime)

      gam = PRIVATE_create_GAM(to_smooth, knots)
      smoothed = PRIVATE_smooth_GAM(gam, pseudotime)

      peak_index = which.max(smoothed)
      peak_pseudotime = max(pseudotime) * (peak_index/100)

      window_start = peak_index-(window_pct/2)
      window_end = peak_index+(window_pct/2)

      window_start = (window_start/100)*max(pseudotime)
      window_end = (window_end/100)*max(pseudotime)
      mean_in = mean(to_smooth[
        to_smooth$pdt>=window_start & to_smooth$pdt <= window_end,
      ]$nc)
      mean_out = mean(to_smooth[
        to_smooth$pdt<window_start | to_smooth$pdt > window_end,
      ]$nc)

      result = data.frame(
        gene=gene,
        peak_pseudotime=peak_pseudotime,
        mean_in_window=mean_in,
        mean_out_window=mean_out,
        ratio=mean_in/mean_out,
        window_start=window_start,
        window_end=window_end,
        deviance_explained = summary(gam)$dev.expl
      )

      return(result)
    },
    BPPARAM = BPPARAM
  )

  return(do.call("rbind", results))
}

#' smooth_gene
#'
#' @description
#' ...
#'
#' @param sce SCE to do the calculations on.
#' @param lineage lineage to obtain TradeSeq smoothers from.
#' @param window_pct the size of the window to consider, as a percentage
#' of the maximum pseudotime value.
#'
#' @return ...
#' @export
#'
#' @concept gene-selection
#' @import mgcv
#'
#' @examples
smooth_gene = function(sce, gene, lineage = 1,
                     window_pct = 10,
                     pseudotime_slot="slingPseudotime_1",
                     knots = 10) {

  if (!(pseudotime_slot %in% colnames(colData(sce)))) {
    stop("Pseudotime slot not in object")
  }

  pseudotime <- colData(sce)[[pseudotime_slot]]
  normalised_counts = normcounts(sce)
  to_smooth = data.frame(nc=normalised_counts[gene,], pdt=pseudotime)

  gam = PRIVATE_create_GAM(to_smooth, knots)
  smoothed = PRIVATE_smooth_GAM(gam, pseudotime)

  return(smoothed)
}


#' plot_gene_peakedness
#'
#' @param sce
#' @param gene_peakedness_df
#' @param gene
#' @param pseudotime_slot
#'
#' @return
#' @export
#'
#' @examples
plot_gene_peakedness = function(sce, gene_peakedness_df, gene, pseudotime_slot="slingPseudotime_1") {

  if (!(pseudotime_slot %in% colnames(colData(sce)))) {
    stop("Pseudotime slot not in object")
  }

  pseudotime <- colData(sce)[[pseudotime_slot]]


  gene_index = which(gene_peakedness_df$gene == gene)

  if (gene_index == 0) {
    stop("Gene not in gene_peakedness_df")
  }

  target = gene_peakedness_df[gene_index,]

  expression = as.data.frame(SingleCellExperiment::normcounts(sce))
  expression = as.data.frame(t(expression))
  expression$pseudotime = c(sce@colData$slingPseudotime_1)
  expression = expression[order(expression$pseudotime),]

  smooth_df = data.frame(
    pdt=((1:100)/100)*max(expression$pseudotime),
    smooth_val=smooth_gene(
      sce = sce,
      gene = gene,
      pseudotime_slot = "slingPseudotime_1",
      knots = 10)
  )
  colnames(smooth_df) = c("pseudotime", target$gene)
  print(target)
  # y = ggplot2::sym(target$gene)
  p = ggplot(NULL, aes_string(x="pseudotime", y= target$gene)) +
    geom_point(data=expression) +
    geom_vline(xintercept=target$peak_pseudotime) +
    geom_vline(xintercept=target$window_start, linetype='dashed') +
    geom_vline(xintercept=target$window_end, linetype='dashed') +
    geom_hline(yintercept=target$mean_in_window, linetype='dashed', color="red") +
    geom_hline(yintercept=target$mean_out_window, linetype='dashed', color="blue") +ggtitle(paste0(target$gene, " (Ratio:", target$ratio, ")")) +
    geom_line(data=smooth_df, color="green")

  return(p)
}

PRIVATE_create_GAM = function(to_smooth, knots) {
  return(
    mgcv::gam(
      nc ~ s(pdt, bs='cr', k=knots),
      data = to_smooth,
      family = "nb"
    )
  )
}

PRIVATE_smooth_GAM = function(gam, pseudotime) {
  return(mgcv::predict.gam(
    gam,
    newdata=data.frame(pdt=((1:100)/100)*max(pseudotime)),
    type="response"
  ))
}

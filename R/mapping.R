#' Map many bulk samples in the same dataframe
#'
#' @concept mapping
#'
#' @param blase_data The [BlaseData] holding the bins.
#' @param bulk_data The whole bulk read matrix as a dataframe. Each row should
#' represent a gene, and each column a sample.
#' @param bootstrap_iterations The number of bootstrapping iterations to run.
#' @param BPPARAM The BiocParallel param for multithreading if desired.
#' Defaults to [BiocParallel::SerialParam()]
#'
#' @import methods
#' @import BiocParallel
#'
#' @return A vector of [MappingResult] objects.
#' @seealso [map_best_bin()]
#' @export
#'
#' @inherit MappingResult-class examples
map_all_best_bins <- function(blase_data, bulk_data,
                              bootstrap_iterations = 200,
                              BPPARAM = BiocParallel::SerialParam()) {

    bulk_data_blase_genes_only = bulk_data[
      rownames(bulk_data) %in% blase_data@genes,
    ]
    dataframes = list()
    for (col_id in colnames(bulk_data)) {
      df = as.data.frame(bulk_data_blase_genes_only[, col_id])
      rownames(df) <- rownames(bulk_data_blase_genes_only)
      colnames(df) <- col_id
      dataframes[[length(dataframes)+1]] = df
    }

    blase_data@pseudobulk_bins = lapply(
      blase_data@pseudobulk_bins,
      function(x) {
        return(x[blase_data@genes,])
      })

    # force execution of lapply
    force(blase_data@pseudobulk_bins)

    results <- BiocParallel::bplapply(
        dataframes,
        function(df) {
            id = colnames(df)[1]
            return(map_best_bin(blase_data, id, df,
                bootstrap_iterations = bootstrap_iterations
            ))
        },
        BPPARAM = BPPARAM
    )

    return(results)
}


#' Map the best matching SC bin for a bulk sample
#'
#' @concept mapping
#'
#' @param blase_data The [BlaseData] holding the bins.
#' @param bulk_id The sample id of the bulk to analyse.
#' @param bulk_data The whole bulk read matrix as a dataframe. Each row should
#' represent a gene, and each column a sample.
#' @param bootstrap_iterations The number of bootstrapping iterations to run.
#'
#' @import methods
#' @import RVAideMemoire
#' @import dqrng
#'
#' @return A [MappingResult] object.
#' @export
#'
#' @inherit MappingResult-class examples
map_best_bin <- function(
    blase_data, bulk_id, bulk_data, bootstrap_iterations = 200) {

    PRIVATE_quality_check_blase_object(blase_data, bulk_data)

    results <- data.frame()
    for (i in blase_data@bins) {
        results <- rbind(results, PRIVATE_map_bin(
            blase_data,
            i,
            bulk_data,
            bulk_id,
            bootstrap_iterations
        ))
    }

    colnames(results) <- c("bin", "correlation", "lower_bound", "upper_bound")

    best_cor <- max(results$correlation)
    best_i <- which.max(results$correlation)

    top2 <- utils::head(
        sort(results$correlation, decreasing = TRUE),
        n = 2
    )
    distance_between_top_2_corrs <- round(top2[1] - top2[2], 4)

    second_best_high_bound <- max(
        results[results$bin != best_i, ]$upper_bound
    )

    best_low_bound <- results[
        results$bin == best_i,
    ]$lower_bound

    confident <- second_best_high_bound < best_low_bound && best_low_bound > 0

    return(methods::new("MappingResult",
        bulk_name = bulk_id,
        best_bin = best_i,
        best_correlation = best_cor,
        top_2_distance = distance_between_top_2_corrs,
        confident_mapping = confident,
        history = results,
        bootstrap_iterations = bootstrap_iterations
    ))
}

PRIVATE_quality_check_blase_object <- function(blase_data, bulk) {
    if (any(is.na(genes(blase_data))) || length(genes(blase_data)) == 0) {
        stop(
            "No genes to map with. ",
            "Please add something to the genes(blase_data) slot."
        )
    }

    if (length(intersect(genes(blase_data), rownames(bulk))) == 0) {
        stop("No genes in genes(blase_data) exist in ",
            "the rows of the bulk dataframe, exiting.")
    }
}

PRIVATE_quality_check_bin <- function(blase_data, i, genes_present) {

    if (ncol(blase_data@pseudobulk_bins[[i]]) <= 1) {
        stop(
            "Not enough cells in bin ",
            as.character(i),
            " to map against, please reduce number of bins (currently ",
            length(blase_data@pseudobulk_bins),
            ") or split by cells"
        )
    }
}

PRIVATE_map_bin <- function(
    blase_data,
    i,
    bulk_data,
    bulk_id,
    bootstrap_iterations) {

    genes_present_in_ref <- blase_data@genes[
        blase_data@genes %in% rownames(blase_data@pseudobulk_bins[[i]])
    ]

    genes_present_in_both <- genes_present_in_ref[
      genes_present_in_ref %in% rownames(bulk_data)
    ]

    if (length(blase_data@genes) != length(genes_present_in_both)) {
      warning(
        "Genes for mapping not all in bulk, using ",
        length(genes_present_in_both),
        " genes available in both reference and bulk."
      )
    }

    if (any(blase_data@genes == bulk_id)) {
      warning("Bulk ID matches a gene, if this fails then check you are",
              "using bulk name and not geneIds:",bulk_id)
    }

    counts_for_top_genes <- bulk_data[genes_present_in_both, as.character(bulk_id)]

    PRIVATE_quality_check_bin(blase_data, i, genes_present_in_both)

    bin_ratios <- blase_data@pseudobulk_bins[[i]][genes_present_in_both, ]

    # all_info_correlation <- stats::cor.test(
    #     unname(Matrix::rowSums(bin_ratios)),
    #     counts_for_top_genes,
    #     method = "spearman",
    #     exact = FALSE
    # )
    # corr_estimate <- unname(all_info_correlation$estimate)
    #
    # bounds <- PRIVATE_bootstrap_bin(
    #     bootstrap_iterations,
    #     bin_ratios,
    #     counts_for_top_genes
    # )
    #
    # return(c(
    #   i,
    #   corr_estimate,
    #   bounds$lower_bound,
    #   bounds$upper_bound
    # ))

    corr = RVAideMemoire::spearman.ci(
      unname(Matrix::rowSums(bin_ratios)),
      counts_for_top_genes,
      nrep = bootstrap_iterations,
      conf.level = 0.95
    )

    return(c(
      i,
      corr$estimate,
      unname(corr$conf.int[1]),
      unname(corr$conf.int[2])
    ))
}

PRIVATE_bootstrap_bin <- function(
    bootstrap_iterations,
    bin_ratios,
    counts_for_top_genes) {
    lower_bound <- 0
    upper_bound <- 1
    if (bootstrap_iterations > 0) {
    #    correlations <- c()
    #    # TODO AM can we speed this up by using apply?
    #     for (j in seq_len(bootstrap_iterations)) {
    #         pseudobulk_sample_names <- sample(
    #             colnames(bin_ratios), ncol(bin_ratios),
    #             replace = TRUE
    #         )
    #         pseudobulk_sample_means <- unname(
    #             Matrix::rowSums(bin_ratios[, pseudobulk_sample_names])
    #         )
    #
    #         corr <- stats::cor.test(
    #             pseudobulk_sample_means,
    #             counts_for_top_genes,
    #             method = "spearman",
    #             exact = FALSE
    #         )
    #         correlations <- c(correlations, unname(corr$estimate))
    #     }
    #
    #     middle_90 <- correlations[{
    #         q <- rank(correlations) / length(correlations)
    #         q < 0.05 | q >= 0.95
    #     }]
    #     lower_bound <- min(middle_90)
    #     upper_bound <- max(middle_90)
    # }

      samples = matrix(colnames(bin_ratios)[
        dqsample.int(
          ncol(bin_ratios),
          ncol(bin_ratios)*bootstrap_iterations,
          replace=TRUE
        )
      ], nrow=bootstrap_iterations)

      speeds = data.frame(getSampleNames=c(), getSampleMeans=c(), corr=c())

      correlations = lapply(seq_len(bootstrap_iterations), function(j) {
        pseudobulk_sample_names = as.vector(samples[j,])

        pseudobulk_sample_means <- unname(
          Matrix::rowSums(bin_ratios[, pseudobulk_sample_names])
        )

        corr <- stats::cor.test(
          pseudobulk_sample_means,
          counts_for_top_genes,
          method = "spearman",
          exact = FALSE
        )

        return(unname(corr$estimate))
      })

        correlations=unlist(correlations, use.names = FALSE)

        middle_95 <- correlations[{
            q <- rank(correlations) / length(correlations)
            q < 0.025 | q >= 0.975
        }]
        lower_bound <- min(middle_95)
        upper_bound <- max(middle_95)
    }

    return(list("lower_bound" = lower_bound, "upper_bound" = upper_bound))
}

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
    results <- BiocParallel::bplapply(
        colnames(bulk_data),
        function(bulk_id) {
            return(map_best_bin(blase_data, bulk_id, bulk_data,
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
    if (any(length(genes_present) != length(blase_data@genes))) {
        warn(
            "Not all genes present in bucket ",
            i,
            " continuing without checking correlation for these genes.\n"
        )
    }

    if (ncol(blase_data@pseudobulk_bins[[i]]) <= 1) {
        stop(
            "Not enough cells in bin ",
            i,
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
    genes_present <- blase_data@genes[
        blase_data@genes %in% rownames(blase_data@pseudobulk_bins[[i]])
    ]
    counts_for_top_genes <- bulk_data[genes_present, as.character(bulk_id)]

    PRIVATE_quality_check_bin(blase_data, i, genes_present)

    bin_ratios <- blase_data@pseudobulk_bins[[i]][genes_present, ]

    all_info_correlation <- stats::cor.test(
        unname(Matrix::rowMeans(bin_ratios)),
        counts_for_top_genes,
        method = "spearman",
        exact = FALSE
    )
    corr_estimate <- unname(all_info_correlation$estimate)

    bounds <- PRIVATE_bootstrap_bin(
        bootstrap_iterations,
        bin_ratios,
        counts_for_top_genes
    )

    return(c(i, corr_estimate, bounds$lower_bound, bounds$upper_bound))
}

PRIVATE_bootstrap_bin <- function(
    bootstrap_iterations,
    bin_ratios,
    counts_for_top_genes) {
    lower_bound <- 0
    upper_bound <- 1
    if (bootstrap_iterations > 0) {
        correlations <- c()
        for (j in seq_len(bootstrap_iterations)) {
            pseudobulk_sample_names <- sample(
                colnames(bin_ratios), ncol(bin_ratios),
                replace = TRUE
            )
            pseudobulk_sample_means <- unname(
                Matrix::rowMeans(bin_ratios[, pseudobulk_sample_names])
            )

            corr <- stats::cor.test(
                pseudobulk_sample_means,
                counts_for_top_genes,
                method = "spearman",
                exact = FALSE
            )
            correlations <- c(correlations, unname(corr$estimate))
        }

        middle_90 <- correlations[{
            q <- rank(correlations) / length(correlations)
            q < 0.05 | q >= 0.95
        }]
        lower_bound <- min(middle_90)
        upper_bound <- max(middle_90)
    }
    return(list("lower_bound" = lower_bound, "upper_bound" = upper_bound))
}

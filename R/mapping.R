#' Map many bulk samples in the same dataframe
#'
#' @concept mapping
#'
#' @param blase_data The [BlaseData] holding the bins.
#' @param bulk_data The whole bulk read matrix as a dataframe. Each row should
#' represent a gene, and each column a sample.
#' @param bootstrap_iterations The number of bootstrapping iterations to run.
#' @param confidence_level The confidence interval to calculate for mappings.
#' Defaults to 90%.
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
                              confidence_level = 0.90,
                              BPPARAM = BiocParallel::SerialParam()) {
    bulk_data_blase_genes_only <- bulk_data[
        rownames(bulk_data) %in% blase_data@genes,
    ]
    dataframes <- list()
    for (col_id in colnames(bulk_data)) {
        df <- as.data.frame(bulk_data_blase_genes_only[, col_id])
        rownames(df) <- rownames(bulk_data_blase_genes_only)
        colnames(df) <- col_id
        dataframes[[length(dataframes) + 1]] <- df
    }

    blase_data@pseudobulk_bins <- lapply(
        blase_data@pseudobulk_bins,
        function(x) {
            return(x[blase_data@genes, ])
        }
    )

    # force execution of lapply
    force(blase_data@pseudobulk_bins)

    results <- BiocParallel::bplapply(
        dataframes,
        function(df) {
            id <- colnames(df)[1]
            return(map_best_bin(blase_data, id, df,
                bootstrap_iterations = bootstrap_iterations,
                confidence_level = confidence_level
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
#' @param confidence_level The confidence interval to calculate for mappings.
#' Defaults to 90%.
#'
#' @import methods
#'
#' @return A [MappingResult] object.
#' @export
#'
#' @inherit MappingResult-class examples
map_best_bin <- function(
    blase_data,
    bulk_id,
    bulk_data,
    bootstrap_iterations = 200,
    confidence_level = 0.90) {
    PRIVATE_quality_check_blase_object(blase_data, bulk_data)

    results <- data.frame()
    for (i in blase_data@bins) {
        results <- rbind(results, PRIVATE_map_bin(
            blase_data,
            i,
            bulk_data,
            bulk_id,
            bootstrap_iterations,
            confidence_level
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
    second_best_high_bound <- max(results[results$bin != best_i, ]$upper_bound)
    best_low_bound <- results[results$bin == best_i, ]$lower_bound
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

#' @keywords internal
PRIVATE_quality_check_blase_object <- function(blase_data, bulk) {
    if (any(is.na(genes(blase_data))) || length(genes(blase_data)) == 0) {
        stop(
            "No genes to map with. ",
            "Please add something to the genes(blase_data) slot."
        )
    }

    if (length(intersect(genes(blase_data), rownames(bulk))) == 0) {
        stop(
            "No genes in genes(blase_data) exist in ",
            "the rows of the bulk dataframe, exiting."
        )
    }
}

#' @keywords internal
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

#' @keywords internal
PRIVATE_map_bin <- function(
    blase_data,
    i,
    bulk_data,
    bulk_id,
    bootstrap_iterations,
    confidence_level) {
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
        warning(
            "Bulk ID matches a gene, if this fails then check you are",
            "using bulk name and not geneIds:", bulk_id
        )
    }

    counts_for_top_genes <- bulk_data[
        genes_present_in_both, as.character(bulk_id)
    ]

    PRIVATE_quality_check_bin(blase_data, i, genes_present_in_both)

    bin_ratios <- blase_data@pseudobulk_bins[[i]][genes_present_in_both, ]

    corr <- PRIVATE_spearman.ci(
        unname(Matrix::rowSums(bin_ratios)),
        counts_for_top_genes,
        nrep = bootstrap_iterations,
        conf.level = confidence_level
    )

    return(c(
        i, corr$estimate, unname(corr$conf.int[1]), unname(corr$conf.int[2])
    ))
}

#' Confidence interval of a Spearman's rank correlation coefficient
#'
#' Computes the confidence interval of a Spearman's rank correlation
#' coefficient by bootstraping. Originally implemented in RVAidemoire
#' Version 0.9-83-7.
#'
#' @param var1 numeric vector (first variable).
#' @param var2 nuermic verctor (second variable).
#' @param nrep number of replicates for bootstrapping.
#' @param conf.level confidence level of the interval.
#'
#' @returns description method name of the test.
#' @returns data.name a character string giving the name(s) of the data.
#' @returns conf.level confidence level.
#' @returns rep number of replicates.
#' @returns estimate Spearman's rank correlation coefficient.
#' @returns conf.int confidence interval.
#'
#' @import boot
#' @keywords internal
#'
PRIVATE_spearman.ci <-
    function(var1, var2, nrep = 1000, conf.level = 0.95) {
        if (length(var1) != length(var2)) {
            stop("'", deparse(substitute(var1)), "' and '",
                deparse(substitute(var2)), "' lengths differ", sep = "")
        }
        data.name <- paste(deparse(substitute(var1)), " and ",
            deparse(substitute(var2)), "\n", nrep, " replicates", sep = "")
        nul <- as.numeric(
            row.names(table(c(which(is.na(var1)), which(is.na(var2)))))
        )
        var1.2 <- if (length(nul) > 0) {
            var1[-nul]
        } else {
            var1
        }
        var2.2 <- if (length(nul) > 0) {
            var2[-nul]
        } else {
            var2
        }
        cor.fun <- function(data, ind) {
            as.numeric(suppressWarnings(
                stats::cor.test(
                  data[ind, 1],
                  data[ind, 2],
                  method = "spearman"
                )$estimate
            ))
        }
        simul <- boot::boot(data.frame(var1.2, var2.2), cor.fun, R = nrep)
        interval <- PRIVATE_.ci(simul$t, conf.level = conf.level)
        attr(interval, "conf.level") <- conf.level
        coeff <- as.numeric(suppressWarnings(
            stats::cor.test(var1, var2, method = "spearman")$estimate
        ))
        names(coeff) <- "rho"
        result <- list(
            method = "Spearman's rank correlation", conf.level = conf.level,
            rep = nrep, data.name = data.name, estimate = coeff,
            conf.int = interval)
        class(result) <- "htest"
        return(result)
    }

#' .ci
#'
#' Originally implemented in RVAidemoire
#' Version 0.9-83-7.
#' @param x data to calculate ci for
#' @param conf.level confidence level to calculate
#'
#' @keywords internal
#' @returns confidence interval results
PRIVATE_.ci <- function(x, conf.level = 0.95) {
    tri <- sort(stats::na.omit(x))
    if (any(!is.finite(tri))) {
        tri <- tri[-which(!is.finite(tri))]
    }
    repet <- length(tri)
    int <- (1 - conf.level) / 2
    if (repet * int < 1) {
        int.inf <- ceiling(repet * int)
    } else {
        int.inf <- floor(repet * int)
    }
    int.sup <- ceiling(repet * (1 - int))
    result <- c("Inf" = tri[int.inf], "Sup" = tri[int.sup])
    return(result)
}

#' Map many bulk samples in the same dataframe
#'
#' @concept mapping
#'
#' @param blase_data The [BlaseData] holding the bins and pseudobulks.
#' @param bulk_data Dataframe. The whole bulk read matrix as a dataframe.
#' Each row should represent a gene, and each column a sample.
#' @param bootstrap_iterations Integer. The number of bootstrapping
#' iterations to run.
#' @param confidence_level Decimal between 0-1. The confidence interval
#' to calculate for mappings. Defaults to 0.9, or 90%.
#' @param BPPARAM The [BiocParallel::BiocParallelParam] for
#' multithreading if desired. Defaults to [BiocParallel::SerialParam()]
#' @param metric Character. The metric to use to compare mappings. One of:
#' 'spearman', 'pearson', 'kendall', 'cosine_similarity', 'euclidean',
#' 'manhattan.'
#'
#' @importFrom methods new
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocParallel bplapply
#'
#' @return A vector of [MappingResult] objects.
#' @seealso [map_best_bin()]
#' @export
#'
#' @inherit MappingResult examples
map_all_best_bins <- function(blase_data, bulk_data,
                              bootstrap_iterations = 200,
                              confidence_level = 0.90,
                              BPPARAM = BiocParallel::SerialParam(),
                              metric = "spearman") {
    bulk_data_blase_genes_only <- bulk_data[
        rownames(bulk_data) %in% genes(blase_data),
    ]
    dataframes <- list()
    for (col_id in colnames(bulk_data)) {
        df <- as.data.frame(bulk_data_blase_genes_only[, col_id])
        rownames(df) <- rownames(bulk_data_blase_genes_only)
        colnames(df) <- col_id
        dataframes[[length(dataframes) + 1]] <- df
    }

    pseudobulk_bins(blase_data) <- lapply(
      pseudobulk_bins(blase_data),
        function(x) {
            return(x[genes(blase_data), ])
        }
    )

    # force execution of lapply
    force(pseudobulk_bins(blase_data))

    results <- BiocParallel::bplapply(
        dataframes,
        function(df) {
            id <- colnames(df)[1]
            return(map_best_bin(blase_data, id, df,
                bootstrap_iterations = bootstrap_iterations,
                confidence_level = confidence_level,
                metric = metric
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
#' @param bulk_id String. The sample id of the bulk to analyse.
#' @param bulk_data Dataframe. The whole bulk read matrix as a dataframe.
#' Each row should represent a gene, and each column a sample.
#' @param bootstrap_iterations Integer. The number of bootstrapping iterations
#' to run.
#' @param confidence_level Decimal between 0-1. The confidence interval to
#' calculate for mappings. Defaults to 90%.
#' @param metric Character. The metric to use to compare mappings. One of:
#' 'spearman', 'pearson', 'kendall', 'cosine_similarity', 'euclidean',
#' 'manhattan.'
#' @param log_data Boolean. When true, bulk and bin values are log2 transformed
#'
#' @return A [MappingResult] object.
#' @export
#'
#' @inherit MappingResult examples
map_best_bin <- function(
    blase_data,
    bulk_id,
    bulk_data,
    bootstrap_iterations = 200,
    confidence_level = 0.90,
    metric = "spearman",
    log_data = FALSE) {

    PRIVATE_quality_check_blase_object(blase_data, bulk_data, bulk_id)

    results <- data.frame()
    for (i in bins(blase_data)) {
        results <- rbind(results, PRIVATE_map_bin(
            blase_data,
            i,
            bulk_data,
            bulk_id,
            bootstrap_iterations,
            confidence_level,
            metric,
            log_data
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
    strong_call <- second_best_high_bound < best_low_bound && best_low_bound > 0

    return(MappingResult(
        bulk_name = bulk_id,
        best_bin = best_i,
        best_correlation = best_cor,
        top_2_distance = distance_between_top_2_corrs,
        strong_mapping = strong_call,
        history = results,
        bootstrap_iterations = bootstrap_iterations,
        metric = metric
    ))
}

#' @keywords internal
PRIVATE_quality_check_blase_object <- function(blase_data, bulk, bulk_id) {

    if (any(genes(blase_data) == bulk_id)) {
      warning(
        "Bulk ID matches a gene, if this fails then check you are",
        "using bulk name and not geneIds:", bulk_id
      )
    }

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
    if (ncol(pseudobulk_bins(blase_data)[[i]]) <= 1) {
        stop(
            "Not enough cells in bin ",
            as.character(i),
            " to map against, please reduce number of bins (currently ",
            length(pseudobulk_bins(blase_data)),
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
    confidence_level,
    metric,
    log_data) {
    genes_in_ref <- genes(blase_data)[
        genes(blase_data) %in% rownames(pseudobulk_bins(blase_data)[[i]])]
    genes_in_both <- genes_in_ref[genes_in_ref %in% rownames(bulk_data)]

    if (length(genes(blase_data)) != length(genes_in_both)) {
        warning("Genes for mapping not all in bulk, using ",
            length(genes_in_both),
            " genes available in both reference and bulk.")
    }
    PRIVATE_quality_check_bin(blase_data, i, genes_in_both)
    counts_for_top_genes <- bulk_data[genes_in_both, as.character(bulk_id)]
    bin_ratios <- pseudobulk_bins(blase_data)[[i]][genes_in_both, ]
    var1 <- unname(Matrix::rowSums(bin_ratios))
    var2 <- counts_for_top_genes
    if (log_data) {
      var1 <- log(var1, base=2)
      var2 <- log(var2, base=2)
    }
    if (metric %in% c("spearman", "pearson", "kendall")) {
      res <- PRIVATE_correlation.ci(
          var1, var2, nrep = bootstrap_iterations,
          conf.level = confidence_level, metric=metric)
    } else if (metric == "cosine_similarity") {
      res <- PRIVATE_cosine_similarity.ci(
        var1, var2, nrep = bootstrap_iterations, conf.level = confidence_level)
    } else if (metric %in% c("euclidean", "manhattan")) {
      res <- PRIVATE_distance.ci(
        var1, var2, nrep = bootstrap_iterations, conf.level = confidence_level,
        metric=metric)
    } else {
      stop("Metric not recognised. Please see documentation for valid options.")
    }
    return(c(
        i, res$estimate, unname(res$conf.int[1]), unname(res$conf.int[2])))
}

#' Confidence interval of a correlation coefficient
#'
#' Computes the confidence interval of a correlation
#' coefficient by bootstraping. Adapted from the implementation of
#' spearman.ci in RVAidemoire
#' Version 0.9-83-7.
#'
#' @param var1 numeric vector (first variable).
#' @param var2 nuermic verctor (second variable).
#' @param nrep number of replicates for bootstrapping.
#' @param conf.level confidence level of the interval.
#' @param metric character from "pearson", "spearman", "kendall", the
#' correlation method to use.
#'
#' @returns description method name of the test.
#' @returns data.name a character string giving the name(s) of the data.
#' @returns conf.level confidence level.
#' @returns rep number of replicates.
#' @returns estimate of correlation coefficient.
#' @returns conf.int confidence interval.
#'
#' @importFrom boot boot
#' @keywords internal
#'
PRIVATE_correlation.ci <-
    function(var1, var2, nrep = 1000, conf.level = 0.95, metric="spearman") {
        if (length(var1) != length(var2)) {
            stop("'", deparse(substitute(var1)), "' and '",
                deparse(substitute(var2)), "' lengths differ",
                sep = "")
        }
        data.name <- paste(deparse(substitute(var1)), " and ",
            deparse(substitute(var2)), "\n", nrep, " replicates", sep = "")
        nul <- as.numeric(
            row.names(table(c(which(is.na(var1)), which(is.na(var2))))))
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
                    method = metric
                )$estimate))
        }
        simul <- boot::boot(data.frame(var1.2, var2.2), cor.fun, R = nrep)
        interval <- PRIVATE_.ci(simul$t, conf.level = conf.level)
        attr(interval, "conf.level") <- conf.level
        coeff <- as.numeric(suppressWarnings(
            stats::cor.test(var1, var2, method = metric)$estimate))

        if (metric == "spearman") {
          names(coeff) <- "rho"
        } else if (metric == "pearson") {
          names(coeff) <- "rho"
        } else if (metric == "kendall") {
          names(coeff) <- "tau"
        }

        result <- list(
            method = metric, conf.level = conf.level, rep = nrep,
            data.name = data.name, estimate = coeff, conf.int = interval)
        class(result) <- "htest"
        return(result)
    }

#' Confidence interval of a distance metric
#'
#' Computes the confidence interval of a Spearman's rank correlation
#' coefficient by bootstraping. Adapted from the implementation
#' of spearman.ci in RVAidemoire
#' Version 0.9-83-7.
#'
#' @param var1 numeric vector (first variable).
#' @param var2 nuermic verctor (second variable).
#' @param nrep number of replicates for bootstrapping.
#' @param conf.level confidence level of the interval.
#' @param metric character from "euclidean", "manhattan", the
#' distance method to use.
#'
#' @returns description method name of the test.
#' @returns data.name a character string giving the name(s) of the data.
#' @returns conf.level confidence level.
#' @returns rep number of replicates.
#' @returns estimate calculated distance (as a negative, so that there is
#' consistency with other methods where a higher value indicates more
#' similarity).
#' @returns conf.int confidence interval.
#'
#' @importFrom boot boot
#' @importFrom stats dist
#' @keywords internal
#'
PRIVATE_distance.ci <-
  function(var1, var2, nrep = 1000, conf.level = 0.95, metric="euclidean") {
    if (length(var1) != length(var2)) {
      stop("'", deparse(substitute(var1)), "' and '",
           deparse(substitute(var2)), "' lengths differ",
           sep = "")
    }
    data.name <- paste(deparse(substitute(var1)), " and ",
                       deparse(substitute(var2)), "\n", nrep, " replicates",
                       sep = ""
    )
    nul <- as.numeric(
      row.names(table(c(which(is.na(var1)), which(is.na(var2))))))

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
    dist.fun <- function(data, ind) {
      as.numeric(suppressWarnings(
        stats::dist(
          matrix(c(data[ind, 1], data[ind, 2]), nrow=2), method = metric)[1]))
    }
    simul <- boot::boot(data.frame(var1.2, var2.2), dist.fun, R = nrep)
    interval <- PRIVATE_.ci(simul$t, conf.level = conf.level)
    attr(interval, "conf.level") <- conf.level
    coeff <- as.numeric(suppressWarnings(
      dist(matrix(c(var1, var2), nrow=2), method = metric)[1]))

    names(coeff) <- "similarity"

    coeff <- coeff*-1
    interval <- interval*-1

    result <- list(
      method = metric, conf.level = conf.level, rep = nrep,
      data.name = data.name, estimate = coeff, conf.int = interval)
    class(result) <- "htest"
    return(result)
  }

#' Confidence interval of cosine similarity
#'
#' Computes the confidence interval of a cosine similarity
#' coefficient by bootstraping. Adapted from the implementation
#' of spearman.ci in RVAidemoire
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
#' @importFrom boot boot
#' @importFrom lsa cosine
#' @keywords internal
#'
PRIVATE_cosine_similarity.ci <-
  function(var1, var2, nrep = 1000, conf.level = 0.95) {

    if (length(var1) != length(var2)) {
      stop("'", deparse(substitute(var1)), "' and '",
           deparse(substitute(var2)), "' lengths differ",
           sep = ""
      )
    }
    data.name <- paste(deparse(substitute(var1)), " and ",
                       deparse(substitute(var2)), "\n", nrep, " replicates",
                       sep = ""
    )
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
    dist.fun <- function(data, ind) {
      as.numeric(suppressWarnings(
        lsa::cosine(data[ind, 1], data[ind, 2])
      ))
    }
    simul <- boot::boot(data.frame(var1.2, var2.2), dist.fun, R = nrep)
    interval <- PRIVATE_.ci(simul$t, conf.level = conf.level)
    attr(interval, "conf.level") <- conf.level
    coeff <- as.numeric(suppressWarnings(
      lsa::cosine(var1, var2)
    ))

    names(coeff) <- "similarity"


    result <- list(
      method = "cosine_similarity", conf.level = conf.level,
      rep = nrep, data.name = data.name, estimate = coeff,
      conf.int = interval
    )
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

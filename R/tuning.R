#' Evaluate n_bins and n_genes for bin mapping
#'
#' @description Will use the n_bins and n_genes implied by the `sce` and
#' `pseudotime_bins_top_n_genes_df` parameters and return quality metrics and
#' an optional chart.
#'
#' @concept tuning
#'
#' @param blase_data The [BlaseData] object to use.
#' @param bootstrap_iterations Iterations for bootstrapping when calculating
#' confident mappings.
#' @param BPPARAM The BiocParallel configuration. Defaults to SerialParam.
#' @param make_plot Whether or not to render the plot showing the correlations
#' for each pseudobulk bin when we try to map the given bin.
#' @param plot_columns How many columns to use in the plot.
#'
#' @return A vector of length 3:
#' * "worst top 2 distance" containing the lowest difference between the
#'  absolute values of the top 2 most correlated bins for each bin.
#'  Higher is better for differentiating.
#' * "mean top 2 distance" containing the mean top 2 distance across the
#'  entire set of genes and bins. Higher is better for differentiation,
#'  but it should matter less than the worst value.
#' * "confident_mapping_pct" - The percent of mappings for this setup which
#'  were annotated as confident by BLASE
#' @export
#'
#' @examples
#' ncells <- 70
#' ngenes <- 100
#' counts_matrix <- matrix(
#'     c(seq_len(3500) / 10, seq_len(3500) / 5),
#'     ncol = ncells,
#'     nrow = ngenes
#' )
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(
#'     normcounts = counts_matrix, logcounts = log(counts_matrix)
#' ))
#' colnames(sce) <- seq_len(ncells)
#' rownames(sce) <- as.character(seq_len(ngenes))
#' sce$cell_type <- c(
#'     rep("celltype_1", ncells / 2),
#'     rep("celltype_2", ncells / 2)
#' )
#'
#' sce$pseudotime <- seq_len(ncells)
#' genelist <- as.character(seq_len(ngenes))
#'
#' # Evaluating created BlaseData
#' blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 10)
#' genes(blase_data) <- genelist[1:20]
#'
#' # Check convexity of parameters
#' evaluate_parameters(blase_data, make_plot = TRUE)
evaluate_parameters <- function(
    blase_data,
    bootstrap_iterations = 200,
    BPPARAM = BiocParallel::SerialParam(),
    make_plot = FALSE,
    plot_columns = 4) {
    results.best_bin <- c()
    results.best_corr <- c()
    results.history <- c()
    results.convexity <- c()
    results.confident_mapping <- c()

    # TODO AM This just pseudobulks every cell
    # pseudobulked_bins <- data.frame(lapply(blase_data@bins, function(i) {
    #     return(Matrix::rowSums(blase_data@pseudobulk_bins[[i]]))
    # }))

    # TODO AM This randomly selects 50% of cells for use on each side
    pseudobulked_bins = NULL
    for (i in blase_data@bins) {
      x <- blase_data@pseudobulk_bins[[i]]
      split = round(runif(ncol(x), 0, 1))
      test = as.matrix(x[,split==1])
      train = as.matrix(x[,split==0])
      test_pseudobulk = Matrix::rowSums(test)
      pseudobulked_bins = cbind(pseudobulked_bins, test_pseudobulk)
      blase_data@pseudobulk_bins[[i]] = train
    }

    colnames(pseudobulked_bins) <- blase_data@bins

    results = map_all_best_bins(blase_data = blase_data,
                                bulk_data = pseudobulked_bins,
                                bootstrap_iterations = bootstrap_iterations,
                                BPPARAM = BPPARAM)

    for (res in results) {
        results.best_bin <- append(results.best_bin, c(res@best_bin))
        results.best_corr <- append(results.best_corr, c(res@best_correlation))
        results.convexity <- append(
            results.convexity, c(res@top_2_distance)
        )
        results.history <- append(results.history, c(res@history))
        results.confident_mapping <- append(results.confident_mapping, c(res@confident_mapping))
    }

    minimum_convexity <- min(results.convexity)
    mean_convexity <- mean(results.convexity)

    # TRUE evaluated as 1
    confident_mapping_pct <- (sum(results.confident_mapping) / length(blase_data@bins)) * 100

    if (make_plot == TRUE) {
        PRIVATE_evaluate_parameters_plots(
            blase_data,
            blase_data@bins,
            results.best_bin,
            results.best_corr,
            results.history,
            results.convexity,
            plot_columns,
            minimum_convexity
        )
    }

    return(c(minimum_convexity, mean_convexity, confident_mapping_pct))
}

PRIVATE_evaluate_parameters_plots <- function(
    blase_data,
    bin_ids,
    results.best_bin,
    results.best_corr,
    results.history,
    results.convexity,
    plot_columns,
    minimum_convexity) {
    plots <- list()

    for (i in seq_len(length(bin_ids))) {
        plots[[i]] <- PRIVATE_plot_history(
            i,
            results.best_bin,
            results.best_corr,
            results.history,
            results.convexity
        )
    }

    gridExtra::grid.arrange(
        top = grid::textGrob(paste(
            length(blase_data@genes),
            "genes and worst convexity:",
            signif(minimum_convexity, 2)
        ), gp = grid::gpar(fontsize = 20, font = 3)),
        grobs = plots,
        ncol = plot_columns
    )
}


#' Identify the Best Parameters For Your Dataset
#'
#' @concept tuning
#'
#' @param x The object to create `BlaseData`` from
#' @param genelist The list of genes to use (ordered by descending goodness)
#' @param bins_count_range The n_bins list to try out
#' @param gene_count_range The n_genes list to try out
#' @param bootstrap_iterations Iterations for bootstrapping when calculating
#' confident mappings.
#' @param BPPARAM The BiocParallel configuration. Defaults to SerialParam.
#' @param verbose Whether to print the n_gene/n_bin combination in progress.
#' Defaults to False.
#' @param ... params to be passed to child functions, see [as.BlaseData()]
#'
#' @return A dataframe of the results.
#' * bin_count: The bin count for this attempt
#' * gene_count: The top n genes to use for this attempt
#' * minimum_convexity: The worst convexity for these parameters
#' * mean_convexity: The mean convexity for these parameters
#' * confident_mapping_pct: The percent of bins which were confidently mapped
#'   to themselves for these parameters. If this value is low, then it is
#'   likely that in real use, few or no results will be confidently mapped.
#'
#' @seealso [plot_find_best_params_results()] for plotting the
#' results of this function.
#'
#' @export
#'
#' @examples
#' ncells <- 70
#' ngenes <- 100
#' counts_matrix <- matrix(
#'     c(seq_len(3500) / 10, seq_len(3500) / 5),
#'     ncol = ncells,
#'     nrow = ngenes
#' )
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(
#'     normcounts = counts_matrix, logcounts = log(counts_matrix)
#' ))
#' colnames(sce) <- seq_len(ncells)
#' rownames(sce) <- as.character(seq_len(ngenes))
#' sce$cell_type <- c(
#'     rep("celltype_1", ncells / 2),
#'     rep("celltype_2", ncells / 2)
#' )
#'
#' sce$pseudotime <- seq_len(ncells)
#' genelist <- as.character(seq_len(ngenes))
#'
#' # Finding the best params for the BlaseData
#' best_params <- find_best_params(
#'     sce, genelist,
#'     bins_count_range = c(10, 20),
#'     gene_count_range = c(20, 50),
#'     pseudotime_slot = "pseudotime",
#'     split_by = "pseudotime_range"
#' )
#' best_params
#' plot_find_best_params_results(best_params)
find_best_params <- function(
    x,
    genelist,
    bins_count_range = c(5, 10, 20, 40),
    gene_count_range = c(10, 20, 40, 80),
    bootstrap_iterations = 200,
    BPPARAM = BiocParallel::SerialParam(),
    verbose = FALSE,
    ...) {
    if (length(genelist) < max(gene_count_range)) {
        stop(
            "Not enough genes provided to meet tuning requests. Provided=",
            length(genelist),
            " wanted=",
            max(gene_count_range)
        )
    }

    results <- data.frame(
        gene_count = c(),
        bin_count = c(),
        minimum_convexity = c(),
        mean_convexity = c(),
        confident_mapping_pct = c()
    )

    for (bin_count in bins_count_range) {
        blase_data <- as.BlaseData(x = x, n_bins = bin_count, ...)

        for (genes_count in gene_count_range) {
            blase_data@genes <- genelist[seq_len(genes_count)]

            if (verbose) {
                message("Bins=", bin_count, " genes=", genes_count)
            }

            res <- evaluate_parameters(
              blase_data,
              bootstrap_iterations,
              BPPARAM,
              make_plot = FALSE
            )
            results <- rbind(
                results,
                data.frame(
                    bin_count = c(bin_count),
                    gene_count = c(genes_count),
                    minimum_convexity = c(res[1]),
                    mean_convexity = c(res[2]),
                    confident_mapping_pct = c(res[3])
                )
            )
        }
    }

    return(results)
}

#' Plot the results of the search for good parameters
#'
#' @concept tuning
#'
#' @param find_best_params_results Results dataframe from [find_best_params()]
#' @param bin_count_colors Optional, custom bin count color scheme.
#' @param gene_count_colors Optional, custom gene count color scheme.
#'
#' @returns A plot showing how convexity changes as n_bins and n_genes
#' are changed. See [find_best_params()] for details on how to interpret.
#'
#' @seealso [find_best_params()]
#'
#' @import viridis
#'
#' @export
#'
#' @inherit find_best_params examples
plot_find_best_params_results <- function(
    find_best_params_results,
    bin_count_colors = viridis::scale_color_viridis(option = "viridis"),
    gene_count_colors = viridis::scale_color_viridis(option = "magma")) {
    gene_count <- ggplot2::sym("gene_count")
    bin_count <- ggplot2::sym("bin_count")
    minimum_convexity <- ggplot2::sym("minimum_convexity")
    mean_convexity <- ggplot2::sym("mean_convexity")
    confident_mapping_pct <- ggplot2::sym("confident_mapping_pct")

    return(gridExtra::grid.arrange(
        # Worst convexity
        ggplot2::ggplot(find_best_params_results, ggplot2::aes(
            x = {{ gene_count }},
            y = {{ minimum_convexity }},
            color = {{ bin_count }}
        )) +
            ggplot2::geom_point() +
            bin_count_colors,
        ggplot2::ggplot(find_best_params_results, ggplot2::aes(
            x = {{ bin_count }},
            y = {{ minimum_convexity }},
            color = {{ gene_count }}
        )) +
            ggplot2::geom_point() +
            gene_count_colors,
        # Mean convexity
        ggplot2::ggplot(find_best_params_results, ggplot2::aes(
            x = {{ gene_count }},
            y = {{ mean_convexity }},
            color = {{ bin_count }}
        )) +
            ggplot2::geom_point() +
            bin_count_colors,
        ggplot2::ggplot(find_best_params_results, ggplot2::aes(
            x = {{ bin_count }},
            y = {{ mean_convexity }},
            color = {{ gene_count }}
        )) +
            ggplot2::geom_point() +
            gene_count_colors,
        # Confident mappings pct
        ggplot2::ggplot(find_best_params_results, ggplot2::aes(
          x = {{ gene_count }},
          y = {{ confident_mapping_pct }},
          color = {{ bin_count }}
        )) +
          ggplot2::geom_point() +
          bin_count_colors,
        ggplot2::ggplot(find_best_params_results, ggplot2::aes(
          x = {{ bin_count }},
          y = {{ confident_mapping_pct }},
          color = {{ gene_count }}
        )) +
          ggplot2::geom_point() +
          gene_count_colors,
        ncol = 2
    ))
}

#' Evaluate Top Genes
#'
#' Shows plots over bins of expression of the top n genes.
#' This is designed to help identify if you have selected
#' genes that vary over the pseudotime you have chosen
#' bins to exist over. Uses the normcounts of the SCE.
#'
#' @concept tuning
#'
#' @param blase_data The [BlaseData] to get bins and expression from.
#' @param n_genes_to_plot The number of genes to plot.
#' @param plot_columns The number of columns to plot the grid with. Best as a
#' divisor of `n_genes_to_plot`.
#'
#' @returns A plot showing the normalised expression of the top genes
#' over pseudotime bins.
#'
#' @export
#'
#' @examples
#' ncells <- 70
#' ngenes <- 100
#' counts_matrix <- matrix(
#'     c(seq_len(3500) / 10, seq_len(3500) / 5),
#'     ncol = ncells,
#'     nrow = ngenes
#' )
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(
#'     normcounts = counts_matrix, logcounts = log(counts_matrix)
#' ))
#' colnames(sce) <- seq_len(ncells)
#' rownames(sce) <- as.character(seq_len(ngenes))
#' sce$cell_type <- c(
#'     rep("celltype_1", ncells / 2),
#'     rep("celltype_2", ncells / 2)
#' )
#'
#' sce$pseudotime <- seq_len(ncells)
#' genelist <- as.character(seq_len(ngenes))
#'
#' # Evaluating created BlaseData
#' blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 10)
#' genes(blase_data) <- genelist[1:20]
#'
#' # Check gene expression over pseudotime
#' evaluate_top_n_genes(blase_data)
evaluate_top_n_genes <- function(
    blase_data,
    n_genes_to_plot = 16,
    plot_columns = 4) {
    if (n_genes_to_plot > length(blase_data@genes)) {
        n_genes_to_plot <- length(blase_data@genes)
    }

    pseudobulked_bins <- data.frame(
        lapply(seq_len(length(blase_data@pseudobulk_bins)), function(i) {
            x <- blase_data@pseudobulk_bins[[i]]
            return(Matrix::rowMeans(x))
        })
    )
    colnames(pseudobulked_bins) <- blase_data@bins

    plots <- list()
    for (i in seq_len(n_genes_to_plot)) {
        plots[[i]] <- PRIVATE_plot_gene_over_bins(
            pseudobulked_bins,
            blase_data@genes[i]
        )
    }

    return(gridExtra::grid.arrange(
        top = grid::textGrob(
            paste(length(blase_data@genes), "genes"),
            gp = grid::gpar(fontsize = 20, font = 3)
        ),
        grobs = plots,
        ncol = plot_columns
    ))
}

PRIVATE_plot_history <- function(i, bin, corr, history, convexity) {
    bin_sym <- ggplot2::sym("bin")
    corr_sym <- ggplot2::sym("correlation")

    return(ggplot2::ggplot(
        as.data.frame(history[(i * 4 - 3):(i * 4)]),
        ggplot2::aes(x = {{ bin_sym }}, y = {{ corr_sym }})
    ) +
        ggplot2::ylim(-1, 1) +
        ggplot2::ggtitle(paste0(
            bin[i],
            " (",
            signif(corr[i], 2),
            ",",
            signif(convexity[i], 2),
            ")"
        )) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = corr[i], linetype = "dashed") +
        ggplot2::geom_vline(xintercept = bin[i], linetype = "dashed"))
}

PRIVATE_plot_gene_over_bins <- function(pseudobulks, gene) {
    bin_sym <- ggplot2::sym("bin")
    expr_sym <- ggplot2::sym("expr")

    expression <- as.data.frame(t(pseudobulks[gene, ]))
    colnames(expression) <- "expr"
    expression$bin <- seq_len(nrow(expression))

    return(ggplot2::ggplot(expression, ggplot2::aes(
        x = {{ bin_sym }},
        y = {{ expr_sym }}
    )) +
        ggplot2::ggtitle(gene) +
        ggplot2::geom_line())
}

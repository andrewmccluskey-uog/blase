#' @title Plot a summary of the mapping result
#'
#' @concept mapping_plots
#'
#' @rdname plot_mapping_result
#' @param x An object to plot on.
#' @param y The [MappingResult] object to plot
#' @param ... additional arguments passed to object-specific methods.
#'
#' @returns A set of plots describing the mapping.
#'
#' @seealso [plot_mapping_result_corr()], [plot_bin_population()]
#'
#' @export
#' @examples
#' counts_matrix <- matrix(
#'     c(seq_len(120) / 10, seq_len(120) / 5),
#'     ncol = 48, nrow = 5
#' )
#' sce <- SingleCellExperiment::SingleCellExperiment(assays = list(
#'     normcounts = counts_matrix, logcounts = log(counts_matrix)
#' ))
#' colnames(sce) <- seq_len(48)
#' rownames(sce) <- as.character(seq_len(5))
#' sce$cell_type <- c(rep("celltype_1", 24), rep("celltype_2", 24))
#'
#' sce$pseudotime <- seq_len(48)
#' blase_data <- as.BlaseData(sce, pseudotime_slot = "pseudotime", n_bins = 4)
#' genes(blase_data) <- as.character(seq_len(5))
#'
#' bulk_counts <- matrix(seq_len(15) * 10, ncol = 3, nrow = 5)
#' colnames(bulk_counts) <- c("A", "B", "C")
#' rownames(bulk_counts) <- as.character(seq_len(5))
#'
#' result <- map_best_bin(blase_data, "B", bulk_counts)
#'
#' # Plot bin
#' sce <- scater::runUMAP(sce)
#' sce <- assign_pseudotime_bins(
#'     sce,
#'     pseudotime_slot = "pseudotime", n_bins = 4
#' )
#' plot_mapping_result(sce, result, group_by_slot = "cell_type")
#'
setGeneric(
    name = "plot_mapping_result",
    signature = c(x = "x", y = "y"),
    def = function(x, y, ...) standardGeneric("plot_mapping_result")
)

#' @rdname plot_mapping_result
#'
#' @param group_by_slot The slot in the
#' [SingleCellExperiment::SingleCellExperiment] to be used as the coloring
#' for the output plot. Passed to [scater::plotUMAP()] as
#' `colour_by`, and will be used to produce a bar chart of
#' populations in the best mapped bin.
#'
#' @import scater
#' @import patchwork
#'
#' @export
setMethod(
    f = "plot_mapping_result",
    signature = c(x = "SingleCellExperiment", y = "MappingResult"),
    definition = function(x, y, group_by_slot) {

      layout = "
      AB
      CD
      EF
      "

      title <- paste0(y@bulk_name, ": Bin ",
                      y@best_bin, ", Cor ", round(y@best_correlation, 4),
                      ", distance ", y@top_2_distance)

      output <- (scater::plotUMAP(x, colour_by = "pseudotime_bin") +
        scater::plotUMAP(x, colour_by = group_by_slot) +
        scater::plotUMAP(
          x[, x$pseudotime_bin == y@best_bin],
          colour_by = "pseudotime_bin"
        ) +
        scater::plotUMAP(
          x[, x$pseudotime_bin == y@best_bin],
          colour_by = group_by_slot
        ) +
        plot_mapping_result_corr(y) +
        plot_bin_population(x, y@best_bin, group_by_slot = group_by_slot) &
        blase_plots_theme()) +
        patchwork::plot_annotation(title = title, theme = blase_titles_theme()) +
        patchwork::plot_layout(design=layout)

      return(output)
    }
)

#' @title Plot the populations of a bin
#'
#' @concept mapping_plots
#'
#' @rdname plot_bin_population
#' @param x An object to plot on.
#' @param bin The bin ID to plot
#' @param ... additional arguments passed to object-specific methods.
#'
#' @returns A ggplot2 object of a plot of population in the given
#' object for this bin.
#'
#' @export
#' @inherit MappingResult-class examples
setGeneric(
    name = "plot_bin_population",
    signature = c(x = "x"),
    def = function(x, bin, ...) standardGeneric("plot_bin_population")
)

#' @rdname plot_bin_population
#'
#' @param group_by_slot The slot in the
#' [SingleCellExperiment::SingleCellExperiment] to be used as the cell
#' type labels.
#'
#' @export
setMethod(
    f = "plot_bin_population",
    signature = c(x = "SingleCellExperiment"),
    definition = function(x, bin, group_by_slot) {
        var1_sym <- ggplot2::sym("Var1")
        freq_sym <- ggplot2::sym("Freq")

        best_bin_population_data <- as.data.frame(
            table(x[, x$pseudotime_bin == bin]@colData[[group_by_slot]])
        )

        return(ggplot2::ggplot(
            best_bin_population_data[best_bin_population_data$Freq > 0, ],
            ggplot2::aes(x = {{ var1_sym }}, y = {{ freq_sym }})
        ) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::ggtitle(paste("Bin", bin)))
    }
)

#' @title Plot a mapping result heatmap
#'
#' @description
#' Plots Spearman's Rho as the fill colour, and adds * if the [MappingResult]
#'  was confidently assigned.
#'
#' @concept mapping_plots
#'
#' @param mapping_result_list A list of [MappingResult] objects to include
#' in the heatmap.
#' @param heatmap_fill_scale The ggplot2 compatible fill gradient scale to
#' apply to the heatmap.
#' @param annotate_confidence Whether to annotate the heatmap with significant results
#' or not, defaults to TRUE.
#' @param annotate_correlation Whether to annotate the heatmap with the
#' correlation of bin to each bulk sample. Defaults to FALSE.
#' or not, defaults to TRUE.
#' @param bin_order The order in which to plot the pseudotime bins along
#' the x-axis.
#'
#' @returns A heatmap showing the correlations of each mapping result across
#' every pseudotime bin.
#'
#' @export
#' @inherit MappingResult-class examples
plot_mapping_result_heatmap <- function(
    mapping_result_list,
    heatmap_fill_scale = ggplot2::scale_fill_gradientn(
        colors = c("blue", "white", "red"), limits = c(-1, 1)
    ),
    annotate_confidence = TRUE,
    annotate_correlation = FALSE,
    bin_order = NULL) {
    if (!all(lapply(mapping_result_list, class) == "MappingResult")) {
        stop("You must provide a list of MappingResult objects only.")
    }

    bulk_results <- data.frame(
        bulk_name = c(),
        pseudotime_bin = c(),
        correlation = c()
    )

    for (mappingResult in mapping_result_list) {
        this_bulk_results <- PRIVATE_get_df_for_this_bulk_to_plot(
            mappingResult, annotate_confidence, annotate_correlation
        )
        bulk_results <- rbind(bulk_results, this_bulk_results)
    }

    bulk_results$bulk_name <- factor(
        bulk_results$bulk_name,
        levels = as.character(unique(bulk_results$bulk_name))
    )

    if (is.null(bin_order)) {
        bin_order <- bulk_results$pseudotime_bin
    }
    bulk_results$pseudotime_bin <- factor(
        bulk_results$pseudotime_bin,
        levels = as.character(unique(bin_order))
    )

    return(PRIVATE_mapping_result_heatmap_plot(
        bulk_results, heatmap_fill_scale,
        annotate_confidence || annotate_correlation
    ))
}

#' @keywords internal
PRIVATE_get_df_for_this_bulk_to_plot <- function(
    mappingResult, annotate_confident, annotate_corr) {
    history <- mappingResult@history


    mapp_corrs = history[, "correlation"]

    confident_mapping = ifelse(
      history[, "bin"] == mappingResult@best_bin &
        rep(
          mappingResult@confident_mapping,
          length(history[, "bin"])
        ),
      "*",
      ""
    )

    labels = rep("", length(history[, "bin"]))
    if (annotate_corr) {
      labels = paste0(labels, round(mapp_corrs,2))
    }
    if (annotate_confident) {
      labels = paste0(labels, confident_mapping)
    }

    return(data.frame(
        bulk_name = rep(mappingResult@bulk_name, nrow(history)),
        pseudotime_bin = history[, "bin"],
        correlation = history[, "correlation"],
        is_best_bin = history[, "bin"] == mappingResult@best_bin,
        label = labels
    ))
}

#' @keywords internal
PRIVATE_mapping_result_heatmap_plot <- function(
    bulk_results, fill_scale, annotate) {
    bulk_name_sym <- ggplot2::sym("bulk_name")
    pseudotime_bin_sym <- ggplot2::sym("pseudotime_bin")
    correlation_sym <- ggplot2::sym("correlation")
    labels_sym <- ggplot2::sym("label")
    is_best_bin_sym <- ggplot2::sym("is_best_bin")

    # Change elow here
    p <- ggplot2::ggplot(bulk_results, ggplot2::aes(
        x = {{ pseudotime_bin_sym }},
        y = {{ bulk_name_sym }},
        fill = {{ correlation_sym }},
        label = {{ labels_sym }},
        color = {{ is_best_bin_sym }}
    )) +
        ggplot2::labs(
            x = "Pseudotime Bin", y = "Bulk Sample", fill = "Correlation"
        ) +
        fill_scale +
        ggplot2::guides(color = "none")

    if (annotate == TRUE) {
        p <- p + ggplot2::geom_tile(ggplot2::aes(
            width = 0.98,
            height = 0.99
        ), linewidth = 0.8) +
            ggplot2::geom_text(fontface = "bold", colour="#000000") +
            ggplot2::scale_color_manual(
                breaks = c(FALSE, TRUE),
                values = c("transparent", "#000000")
            )
    } else {
        p <- p + ggplot2::geom_tile() +
            ggplot2::scale_color_manual(
                breaks = c(FALSE, TRUE),
                values = c("transparent", "transparent")
            )
    }
    return(p)
}

#' @title Plot a mapping result's correlation
#'
#' @description
#' Plots the mapping results correlations with each pseudotime bin
#'
#' @concept mapping_plots
#'
#' @param mapping_result A [MappingResult] object to plot the correlations for.
#'
#' @returns A [ggplot2] object of the plot
#'
#' @export
#' @inherit MappingResult-class examples
plot_mapping_result_corr <- function(mapping_result) {
    bin_sym <- ggplot2::sym("bin")
    correlation_sym <- ggplot2::sym("correlation")
    upper_bound_sym <- ggplot2::sym("upper_bound")
    lower_bound_sym <- ggplot2::sym("lower_bound")
    return(ggplot2::ggplot(mapping_result@history, ggplot2::aes(
        x = {{ bin_sym }}, y = {{ correlation_sym }}
    )) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(
            yintercept = mapping_result@best_correlation,
            linetype = "dashed"
        ) +
        ggplot2::geom_vline(
            xintercept = mapping_result@best_bin,
            linetype = "dashed"
        ) +
        ggplot2::geom_line(
            ggplot2::aes(y = {{ lower_bound_sym }}),
            linetype = "dotted"
        ) +
        ggplot2::geom_line(
            ggplot2::aes(y = {{ upper_bound_sym }}),
            linetype = "dotted"
        ))
}

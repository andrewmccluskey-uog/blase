#' @keywords internal
blase_plots_theme <- function() {
    return(ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12)
    ))
}

#' @keywords internal
blase_titles_theme <- function() {
    return(ggplot2::theme(plot.title = ggplot2::element_text(size = 18)))
}

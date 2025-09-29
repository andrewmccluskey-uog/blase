#' @keywords internal
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
blase_plots_theme <- function() {
    return(ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12)
    ))
}

#' @keywords internal
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
blase_titles <- function() {
    return(ggplot2::theme(plot.title = ggplot2::element_text(size = 18)))
}

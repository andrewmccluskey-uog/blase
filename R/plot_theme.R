#' keywords internal
blase_plots_theme <- function(rotate_x_labels = FALSE) {


  if (rotate_x_labels) {
    return(ggplot2::theme(
      plot.title = element_text(face="bold", size = 12),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
  } else {
    return(ggplot2::theme(
      plot.title = element_text(face="bold", size = 12)))
  }

}

blase_titles_theme <- function() {
  return(ggplot2::theme(plot.title = ggplot2::element_text(size = 18)))
}

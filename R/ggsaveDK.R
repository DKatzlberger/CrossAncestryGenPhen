#' Save ggplot with optional legend removal and sensible defaults
#'
#' @param plot A ggplot object.
#' @param filename File path to save to.
#' @param legend Logical. If FALSE, removes the legend before saving.
#' @param width Width in inches. Default is 7.
#' @param height Height in inches. Default is 7.
#' @param dpi Dots per inch. Default is 300.
#' @param ... Additional arguments passed to ggsave().
#'
#' @return Invisibly returns the saved plot.
#' 
#' @importFrom ggplot2 ggsave theme
#' @export
ggsaveDK <- function(
  filename, 
  plot, 
  legend = TRUE,
  width = 7, 
  height = 7, 
  dpi = 300, 
  ...
) {

  # Remove legend if needed
  if (!legend) {
    plot <- plot + ggplot2::theme(legend.position = "none")
  }
  
  # Save with defaults
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    ...
  )
  
  invisible(plot)
}

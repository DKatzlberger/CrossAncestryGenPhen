#' Save ggplot with optional legend removal and sensible defaults
#'
#' @param plot A ggplot object.
#' @param file File path to save to.
#' @param guides Logical. If FALSE, removes the legend before saving.
#' @param trimmed Logical. If TRUE, removes plot margins before saving (for tight SVG export).
#' @param units Units for width/height. Default 'cm'.
#' @param width Width of saved plot. Default 8.
#' @param height Height of saved plot. Default 8.
#' @param bg Background color for saving (e.g. "white" or "transparent").
#' @param ... Additional arguments passed to [ggplot2::ggsave()].
#'
#' @return Invisibly returns the saved plot.
#'
#' @importFrom ggplot2 ggsave theme
#' @export
ggsaveDK <- function(
  file,
  plot,
  guides = TRUE,
  trimmed = TRUE,
  units = "cm",
  width = 8,
  height = 8,
  bg = "transparent",
  ...
) {

  ## --- Helper: remove guides ---
  noGuides <- function() {
    ggplot2::theme(
      legend.position = "none",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank()
    )
  }

  ## --- Helper: trim margins ---
  trimMargins <- function() {
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
  }

  ## --- Apply optional adjustments ---
  if (trimmed) {
    plot <- plot + trimMargins()
  }

  if (!guides) {
    plot <- plot + noGuides()
    file <- gsub("(\\.[a-z]+$)", "_noGuide\\1", file)
  }

  ## --- Save ---
  ggplot2::ggsave(
    filename = file,
    plot = plot,
    units = units,
    width = width,
    height = height,
    bg = bg,
    ...
  )

  invisible(plot)
}

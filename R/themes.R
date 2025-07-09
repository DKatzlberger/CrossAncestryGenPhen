#' Nature-Inspired Font Sizes Theme (Internal)
#'
#' A ggplot2 theme that sets small, uniform base font sizes for a minimal,
#' nature-journal-style aesthetic. Designed for compact, consistent typography.
#'
#' @param base_size Base font size for all text elements (default: 5).
#' @return A ggplot2 theme object.
#' @keywords internal
#' @noRd
theme_nature_fonts <- function(base_size = 8) {
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = base_size),
    axis.title = ggplot2::element_text(size = base_size),
    plot.title = ggplot2::element_text(size = base_size, hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0.5),
    legend.title = ggplot2::element_text(size = base_size),
    legend.text = ggplot2::element_text(size = base_size),
    strip.text = ggplot2::element_text(size = base_size),
    plot.caption = ggplot2::element_text(size = base_size, hjust = 0)
  )
}

#' Small Legend Theme (Internal)
#'
#' A helper ggplot2 theme that minimizes legend spacing and key size,
#' useful for compact plots or figure panels.
#'
#' @param base_size Base text size to scale legend keys (default: 5).
#' @param ... Additional arguments passed to [ggplot2::theme()].
#'
#' @return A ggplot2 theme object with small legend spacing.
#' @keywords internal
#' @noRd
theme_small_legend <- function(base_size = 8, ...) {
  key_size_pt <- base_size * 1.3  # visually balances with text height

  ggplot2::theme(
    legend.key.spacing = ggplot2::unit(0.2, "pt"),
    legend.key.height = ggplot2::unit(key_size_pt, "pt"),
    legend.key.width = ggplot2::unit(key_size_pt, "pt"),
    legend.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5),
    ...
  )
}


#' @keywords internal
#' @noRd
theme_white_background <- function(...) {
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(color = "black", linewidth = 0.3),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = ggplot2::unit(0.5, "mm"),
    strip.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.text = ggplot2::element_text(color = "black"),
    strip.placement = "inside",
    ...
  )
}


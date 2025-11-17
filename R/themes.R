#' CrossAncestryGenPhenPortability Theme
#'
#' A ggplot2 theme that sets small, uniform base font sizes for a minimal,
#' nature-journal-style aesthetic. Designed for compact, consistent typography.
#'
#' @param base_size Base font size for all text elements (default: 6).
#' @param legend_key Legend key size multiplicator of base_size (default: 1).
#' @param show_axis Logical; whether to show axis lines and ticks.
#' @param show_facets Logical; whether to show facet strip text.
#' @param show_grid Logical; whether to show major grid lines.
#' @param show_boarders Logical; whether to show panel borders.
#' @param rotate Numeric; angle (in degrees) to rotate x-axis labels (e.g. 45, 90).
#' @param base_family Base font family (default: "Arial").
#' @param ... Additional theme parameters.
#'
#' @return A ggplot2 theme object.
#' @export
theme_CrossAncestryGenPhen <- function(
  base_size = 6,
  legend_key = 1.5,
  show_axis = TRUE,
  show_facets = TRUE,
  show_grid = FALSE,
  show_boarders = FALSE,
  rotate = NULL,
  base_family = "Arial",
  ...
) {

  ## --- Unit conversion helpers ---
  mm_to_pt <- function(mm) mm / 0.352778
  pt_to_mm <- function(pt) pt * 0.352778

  ## Base theme
  p <- ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.placement  = "inside",

    axis.line  = ggplot2::element_line(color = "black", linewidth = 0.5),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = ggplot2::unit(1, "mm"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = if (show_grid) ggplot2::element_line(color = "grey80", linewidth = 0.3) else ggplot2::element_blank(),

    axis.text  = ggplot2::element_text(size = base_size, family = base_family),
    axis.title = ggplot2::element_text(size = base_size, family = base_family),
    plot.title = ggplot2::element_text(size = base_size, hjust = 0.5, family = base_family),
    plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0.5, family = base_family),
    strip.text = ggplot2::element_text(size = base_size, family = base_family),

    legend.title = ggplot2::element_text(size = base_size, family = base_family, margin = ggplot2::margin(b = 2, unit = "mm")),
    legend.text = ggplot2::element_text(size = base_size, family = base_family),
    ...
  )

  ## --- Legend sizing ---
  key_size_mm <- pt_to_mm(base_size * legend_key)
  p <- p + ggplot2::theme(
    legend.key.height  = ggplot2::unit(key_size_mm, "mm"),
    legend.key.width   = ggplot2::unit(key_size_mm, "mm"),
    legend.key.spacing = ggplot2::unit(0.5, "mm"),
    legend.margin      = ggplot2::margin(t = 2, b = 2, l = 2, r = 2, unit = "mm"),
    legend.box.spacing = ggplot2::unit(1, "mm"),
    legend.background  = ggplot2::element_rect(fill = NA, color = NA)
  )


  ## --- Adaptive facets ---
  if (!show_facets) p <- p + ggplot2::theme(strip.text = ggplot2::element_blank())

  ## --- Adaptive boarders ---
  if (show_boarders) { p <- p + ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line    = ggplot2::element_blank()
    )
  }

  ## --- Adaptive axis ---
  if (!show_axis) {
    p <- p + ggplot2::theme(
      axis.line   = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(), 
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0.5, b = 0, unit = "mm")), # pull x labels closer
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 0.5, l = 0, unit = "mm"))  # pull y labels closer
    )
  }

  ## --- Adaptive rotation ---
  if (!is.null(rotate)) {
    angle <- rotate %% 360
    # Determine justification heuristically
    if (angle == 0) {
      hjust <- 0.5; vjust <- 0.5
    } else if (angle == 90) {
      hjust <- 1; vjust <- 0.5
    } else if (angle == 270) {
      hjust <- 0; vjust <- 0.5
    } else if (angle > 0 && angle < 90) {
      hjust <- 1; vjust <- 1
    } else if (angle > 90 && angle < 180) {
      hjust <- 1; vjust <- 0
    } else if (angle > 180 && angle < 270) {
      hjust <- 0; vjust <- 0
    } else { # e.g., 315°
      hjust <- 0; vjust <- 1
    }

    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = angle,
        hjust = hjust,
        vjust = vjust,
        family = base_family
      )
    )
  }

  return(p)
}



# #' Small Legend Theme (Internal)
# #'
# #' A helper ggplot2 theme that minimizes legend spacing and key size,
# #' useful for compact plots or figure panels.
# #'
# #' @param base_size Base text size to scale legend keys (default: 5).
# #' @param ... Additional arguments passed to [ggplot2::theme()].
# #' @return A ggplot2 theme object.
# #' @export
# theme_small_legend <- function(base_size = 6, ...) {
#   key_size_pt <- base_size * 1.1 

#   ggplot2::theme(
#     legend.key.spacing = ggplot2::unit(0.2, "pt"),
#     legend.key.height = ggplot2::unit(key_size_pt, "pt"),
#     legend.key.width = ggplot2::unit(key_size_pt, "pt"),
#     legend.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5),
#     legend.title = ggplot2::element_text(margin = ggplot2::margin(b = 5)),
#     ...
#   )
# }


# #' A clean white ggplot2 theme with optional facet labels
# #'
# #' Custom theme with white background, minimal gridlines, and optional control
# #' over facet label display.
# #'
# #' @param show_facets Logical. If `TRUE`, facet strip labels are shown. If `FALSE`, facet strip text is removed.
# #' @param ... Additional arguments passed to [ggplot2::theme()].
# #' @return A ggplot2 theme object.
# #' @export
# theme_white_background <- function(show_facets = TRUE, ...) {
#   ggplot2::theme(
#     panel.background = ggplot2::element_rect(fill = "white", color = NA),
#     plot.background = ggplot2::element_rect(fill = "white", color = NA),
#     panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank(),
#     axis.line = ggplot2::element_line(color = "black", linewidth = 0.3),
#     axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
#     axis.ticks.length = ggplot2::unit(0.5, "mm"),
#     strip.background = ggplot2::element_rect(fill = "white", color = NA),
#     strip.text = if (show_facets) ggplot2::element_text(color = "black") else ggplot2::element_blank(),
#     strip.placement = "inside",
#     ...
#   )
# }


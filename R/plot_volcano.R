#' Volcano Plot
#'
#' Creates a volcano plot with capped -log10(p-values) on the y-axis.
#' Supports effect size threshold lines, feature annotations, and flexible
#' color and shape aesthetics.
#'
#' @param data Data frame with variables to plot.
#' @param x_var String. Column name for x-axis (e.g., effect size).
#' @param y_var String. Column name for y-axis (p-values to transform).
#' @param color_var String. Column name mapped to point colors.
#' @param shape_var String. Column name mapped to point shapes.
#' @param features Optional vector. Feature names to annotate.
#' @param sig_thr Numeric. Significance threshold for horizontal line.
#' @param effect_thr Numeric. Effect size threshold for vertical lines.
#' @param title String. Plot title.
#' @param x_label String. Custom x-axis label.
#' @param y_label String. Custom y-axis label.
#' @param color_label String. Legend title for color.
#' @param shape_label String. Legend title for shape.
#' @param point_size Numeric. Size of points. Default is 1.
#'
#' @return A ggplot2 object.
#'
#' @import ggplot2
#' @import ggrepel
#' @export
plot_volcano <- function(
  data, 
  x_var, 
  y_var, 
  color_var = NULL,
  shape_var = NULL,
  features = NULL,
  sig_thr = NULL, 
  effect_thr = NULL, 
  title = NULL, 
  x_label = NULL,
  y_label = NULL,
  color_label = NULL,
  shape_label = NULL,
  point_size = 1
) {
  
  df <- as.data.frame(data)

  ## --- Dynamic mapping ---
  dyn_mapping <- aes(
    !!!rlang::syms(
      c(
        x = x_var,
        y = y_var,
        if (!is.null(color_var)) c(color = color_var),
        if (!is.null(shape_var)) c(shape = shape_var)
      )
    )
  )

  ## --- Annotate features if provided and 'feature' exists ---
  if (!is.null(features) && "feature" %in% colnames(df)) {
    df$annotate <- ifelse(df$feature %in% features, df$feature, NA)
  } else {
    df$annotate <- NA
  }


  ## --- Base plot ---
  p <- ggplot(
      data = df, 
      mapping = dyn_mapping
    ) +
    geom_point(
      size = point_size
    )
  
  ## --- Axis scale ---
  neglog10_trans <- scales::trans_new(
    name      = "neglog10",
    transform = function(x) -log10(x + 1e-12),
    inverse   = function(x) 10^(-x),
    domain    = c(1e-300, 1),
    breaks    = scales::trans_breaks(function(x) -log10(x), function(x) 10^(-x)),
    format    = scales::scientific_format(digits = 1)
  )

  p <- p + scale_y_continuous(
    trans = neglog10_trans
  )


  ## --- Threshold lines ---
  if (!is.null(sig_thr)) {
    p <- p + geom_hline(
      yintercept = sig_thr,
      linetype = "dashed", 
      linewidth = 0.3,
      color = "blue"
    )
  }

  if (!is.null(effect_thr)) {
    p <- p + geom_vline(
      xintercept = c(
        -effect_thr, 
        effect_thr
      ), 
      linetype = "dashed", 
      linewidth = 0.3,
      color = "blue"
    )
  }

  ## --- Feature annotation ---
  if (!is.null(features) && "feature" %in% colnames(df)) {
    p <- p + ggrepel::geom_text_repel(
      data = subset(df, !is.na(annotate)),
      aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
        label = annotate
      ),
      size = 2,
      min.segment.length = 0,
      segment.size = 0.3,
      max.overlaps = Inf,
      force = 2,
      box.padding = 0.5,
      bg.color = "white",
      bg.r = 0.2,
      inherit.aes = FALSE
    )
  }

  ## --- Final styling ---
  p <- p +
    labs(
      title = title,
      x = ifelse(is.null(x_label), x_var, x_label),
      y = ifelse(is.null(y_label), y_var, y_label),
      color = color_label,
      shape = shape_label
    ) +
    theme_CrossAncestryGenPhen()

  return(p)
}

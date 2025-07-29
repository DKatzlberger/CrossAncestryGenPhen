#' True vs Observed Interaction Effect Plot
#'
#' Creates a scatter plot comparing true and observed interaction effects.
#' Supports threshold band, feature labels, and flexible aesthetics.
#'
#' @param data Data frame with variables for axes and aesthetics.
#' @param x_var String. Column name for x-axis (true effect).
#' @param y_var String. Column name for y-axis (observed effect).
#' @param color_var String. Column name mapped to point colors.
#' @param shape_var String. Column name mapped to point shapes.
#' @param effect_thr Numeric. Adds ribbon band around null effect line.
#' @param features Vector of feature names to annotate.
#' @param title String. Plot title.
#' @param x_label String. Custom x-axis label.
#' @param y_label String. Custom y-axis label.
#' @param color_label String. Legend title for color.
#' @param shape_label String. Legend title for shape.
#' @param fill_label String. Legend title for threshold band.
#' @param caption String. Caption text.
#' @param point_size Numeric. Size of points.
#'
#' @return A ggplot2 object.
#'
#' @import ggplot2
#' @import ggrepel
#' @export
plot_sim_true_vs_obs_interaction_effect <- function(
  data,
  x_var,
  y_var,
  color_var = NULL,
  shape_var = NULL,
  effect_thr = NULL,
  features = NULL,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  color_label = NULL,
  shape_label = NULL,
  fill_label = NULL,
  caption = NULL,
  point_size = 1
) {

  df <- as.data.frame(data)

  # Dynamic mapping
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

  # Annotate features if provided and 'feature' exists
  if (!is.null(features) && "feature" %in% colnames(df)) {
    df$annotate <- ifelse(df$feature %in% features, df$feature, NA)
  } else {
    df$annotate <- NA
  }

  # Abline for 0 interactions
  slope <- 0
  intercept <- 0
  # Span the full x range you want
  x_min <- min(df[[x_var]], na.rm = TRUE)
  x_max <- max(df[[x_var]], na.rm = TRUE)
  x_seq <- seq(x_min, x_max, length.out = 200)
  y_seq <- slope * x_seq + intercept

  line_df <- data.frame(
    x = x_seq, 
    y = y_seq
  )

  # Ribbon for threshold
  if (!is.null(effect_thr)) {
    line_df$ymin <- y_seq - effect_thr
    line_df$ymax <- y_seq + effect_thr
    line_df$Band <- factor(paste("±", effect_thr, "log2FC"))
  }

  # Base plot
  p <- ggplot(
    data = df,
    mapping = dyn_mapping
    ) +
    geom_point(
      size = point_size
    )

    # Ribbon if effect thr is defined (biological interpretation)
    if (!is.null(effect_thr)) {
      p <- p + geom_line(
        data = line_df,
        mapping = aes(
          x = x, 
          y = y,
        ),
        color = "grey50",
        linewidth = 0.3,
        inherit.aes = FALSE
      ) +
      geom_ribbon(
        data = line_df,
        mapping = aes(
          x = x, 
          ymin = ymin, 
          ymax = ymax,
          fill = Band,
          group = Band
        ),
        alpha = 0.2,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(
        values = setNames(
          "grey50", 
          unique(line_df$Band)
        )
      )
    } 

    # Final styling of the plot
    p <- p + labs(
      title = title,
      x = ifelse(is.null(x_label), x_var, x_label),
      y = ifelse(is.null(y_label), y_var, y_label),
      color = color_label,
      shape = shape_label,
      fill = ifelse(is.null(fill_label), "Effect threshold", fill_label),
      caption = caption
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  # Add labels if features are provided
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

  return(p)
}

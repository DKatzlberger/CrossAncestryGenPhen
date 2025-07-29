#' Scatter plot: observed vs. true interaction effects (diff-of-diff)
#'
#' Visualize how observed interaction effects align with true effects in a
#' simulation. Points on the diagonal show perfect estimation; a horizontal
#' ribbon indicates a filtering threshold for the observed effect.
#'
#' @param fX Data frame for ancestry X features.
#' @param fY Data frame for ancestry Y features.
#' @param x_var Column in fX/fY for true ancestry effect.
#' @param y_var Column in fX/fY for observed ancestry effect.
#' @param effect_thr Optional numeric threshold for ribbon band.
#' @param features Optional vector of feature names to annotate.
#' @param title Plot title.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param caption Optional plot caption.
#' @param point_size Size of scatter points.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_hline geom_ribbon
#' @importFrom ggplot2 geom_line scale_color_manual scale_shape_manual
#' @importFrom ggplot2 scale_fill_manual labs theme
#' @importFrom ggrepel geom_text_repel
plot_sim_interaction_effect <- function(
  fX,
  fY,
  x_var,
  y_var,
  effect_thr = NULL,
  features = NULL,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  caption = NULL,
  point_size = 1
) {

  X_a <- unique(fX$ancestry)
  Y_a <- unique(fY$ancestry)

  # Combined effect data frame
  combined_df <- data.frame(
    x = fY[[x_var]] - fX[[x_var]],
    y = fY[[y_var]] - fX[[y_var]],
    X_true_effect = fX$true_log2FC,
    Y_true_effect = fY$true_log2FC,
    is_DE_X = as.logical(fX$is_DE),
    is_DE_Y = as.logical(fY$is_DE),
    feature = fX$feature,
    stringsAsFactors = FALSE
  )

  # Annotate features if provided
  combined_df$annotate <- if (!is.null(features)) {
    ifelse(
      combined_df$feature %in% features, 
      combined_df$feature, 
      NA
    )
  } else {
    NA
  }

  # Fill by true DEGs
  X_status <- paste(X_a, "only")
  Y_status <- paste(Y_a, "only")

  combined_df$DEG_status <- with(
    combined_df, ifelse(
      is_DE_X & is_DE_Y, "Both",
      ifelse(is_DE_X, X_status,
        ifelse(is_DE_Y, Y_status, "None")
      )
    )
  )

  combined_df$DEG_status <- factor(
    combined_df$DEG_status,
    levels = c("Both", X_status, Y_status, "None")
  )

  # Order of point layers
  plot_order <- ifelse(
    combined_df$DEG_status == "Both", 3,
    ifelse(combined_df$DEG_status == "None", 1, 2)
  )
  combined_df <- combined_df[order(plot_order), ]

  # Shape by interactions 
  true_interaction <- combined_df$Y_true_effect - combined_df$X_true_effect
  is_DE <- as.numeric(abs(true_interaction) > 0)
  combined_df$Interaction_status <- factor(
    is_DE,
    levels = c(0, 1),
    labels = c("No interaction", "Interaction")
  )

  # Abline for 0 interactions
  slope <- 0
  intercept <- 0
  # Span the full x range you want
  x_min <- min(combined_df$x, na.rm = TRUE)
  x_max <- max(combined_df$x, na.rm = TRUE)
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
    data = combined_df,
    mapping =  aes(
      x = x, 
      y = y
      )
    ) +
    geom_point(
      mapping = aes(
        color = DEG_status,
        shape = Interaction_status
      ),
      size = point_size
    ) +
    scale_color_manual(
      values = setNames(
        c(
          "purple", 
          "#1f77b4", 
          "#ff7f0e", 
          "gray80"
        ),
        c(
          "Both", 
          X_status, 
          Y_status, 
          "None"
        )
      ),
      breaks = c(
        "Both", 
        X_status, 
        Y_status, 
        "None"
      )
    ) +
    scale_shape_manual(
      values = c(
        "No interaction" = 1, 
        "Interaction" = 17
      ),
      drop = FALSE
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
      caption = caption,
      fill = "Interaction threshold",
      color = "Main effect",
      shape = "Interaction effect"
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  # Add labels if features are provided
  if (!is.null(features)) {
    p <- p + ggrepel::geom_text_repel(
      data = subset(combined_df, !is.na(annotate)),
      aes(
        x = x,
        y = y,
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

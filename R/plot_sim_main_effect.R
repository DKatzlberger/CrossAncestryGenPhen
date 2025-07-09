#' Scatter Plot of Realized log2 Fold-Change X vs Y
#'
#' Create a scatter plot comparing realized log2FC for ancestry X vs. Y.
#'
#' @param fX Data frame for ancestry X features.
#' @param fY Data frame for ancestry Y features.
#' @param x_var Column name in `fX` to use for X-axis effect.
#' @param y_var Column name in `fY` to use for Y-axis effect.
#' @param features Optional character vector of feature names to annotate.
#' @param title Plot title (optional).
#' @param x_label X-axis label (optional).
#' @param y_label Y-axis label (optional).
#' @param point_size Size of points (default: 1).
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggrepel geom_text_repel
plot_sim_main_effect <- function(
  fX,
  fY,
  x_var,
  y_var,
  features = NULL,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  point_size = 1
) {
  stopifnot(identical(fX$feature, fY$feature))

  # Unique intended FCs for X and Y
  X_intended <- unique(fX$intended_log2FC)
  Y_intended <- unique(fY$intended_log2FC)

  X_a <- unique(fX$ancestry)
  Y_a <- unique(fY$ancestry)

  # Combined effect data frame
  combined_df <- data.frame(
    X_effect = fX[[x_var]],
    Y_effect = fY[[y_var]],
    is_DE_X = as.logical(fX$is_DE),
    is_DE_Y = as.logical(fY$is_DE),
    feature = fX$feature,
    stringsAsFactors = FALSE
  )

  # Annotate features if provided
  combined_df$annotate <- if (!is.null(features)) {
    ifelse(combined_df$feature %in% features, combined_df$feature, NA)
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

  # Color by true interactions
  combined_df$Interaction_status <- with(
    combined_df, ifelse(
      xor(is_DE_X, is_DE_Y), "Interaction", "No interaction"
    )
  )

  combined_df$Interaction_status <- factor(
    combined_df$Interaction_status,
    levels = c("Interaction", "No interaction")
  )

  # Lines for intended only
  vline_type <- paste(X_a, "intended effect")
  vline_df <- data.frame(x = c(0, 0), y = c(0, 0), type = vline_type)

  hline_type <- paste(Y_a, "intended effect")
  hline_df <- data.frame(x = c(0, 0), y = c(0, 0), type = hline_type)

  # Plot
  p <- ggplot(
    data = combined_df,
    aes(x = X_effect, y = Y_effect)
  ) +
    geom_point(
      aes(
        color = DEG_status,
        shape = Interaction_status
      ),
      size = point_size
    ) +
    geom_vline(
      xintercept = c(
        -X_intended,
        X_intended
      ),
      linetype = "dashed",
      linewidth = 0.3,
      color = "red",
      show.legend = FALSE
    ) +
    geom_hline(
      yintercept = c(
        -Y_intended,
        Y_intended
      ),
      linetype = "dashed",
      linewidth = 0.3,
      color = "red",
      show.legend = FALSE
    ) +
    geom_line(
      data = vline_df,
      aes(
        x = x, 
        y = y, 
        linetype = type
      ),
      color = "red"
    ) +
    geom_line(
      data = hline_df,
      aes(
        x = x, 
        y = y, 
        linetype = type
      ),
      color = "red"
    ) +
    scale_linetype_manual(
      name = "Intended effect size",
      values = setNames(
        c("dashed", "dashed"),
        c(vline_type, hline_type)
      )
    ) +
    scale_color_manual(
      values = setNames(
        c("purple", "#1f77b4", "#ff7f0e", "gray80"),
        c("Both", X_status, Y_status, "None")
      ),
      breaks = c("Both", X_status, Y_status, "None")
    ) +
    labs(
      title = title,
      x = ifelse(is.null(x_label), paste(X_a, x_var), x_label),
      y = ifelse(is.null(y_label), paste(Y_a, y_var), y_label),
      color = "Intended DEGs",
      shape = "Intended interaction"
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  # Add labels if features are provided
  if (!is.null(features)) {
    p <- p + ggrepel::geom_text_repel(
      data = subset(combined_df, !is.na(annotate)),
      aes(
        x = X_effect,
        y = Y_effect,
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

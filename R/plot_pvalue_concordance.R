#' Plot concordance of -log10 p-values between two methods
#'
#' Creates a scatter plot comparing -log10 transformed p-values from two
#' statistical methods. Optionally colors points by significance and facets
#' by a grouping variable. Correlations and regression lines are shown per
#' facet.
#'
#' @param data Data frame with p-values and optional grouping/significance.
#' @param x_var Character. Column name for x-axis p-values.
#' @param y_var Character. Column name for y-axis p-values.
#' @param x_sig_source Optional. Column for x significance p-values.
#' @param y_sig_source Optional. Column for y significance p-values.
#' @param facet_col Optional. Column to facet by.
#' @param facet_levels Optional. Levels to include from facet_col.
#' @param x_label Character. X-axis label. Defaults to -log10(x_var).
#' @param y_label Character. Y-axis label. Defaults to -log10(y_var).
#' @param title Character. Optional plot title.
#' @param sig_thr Numeric. P-value threshold for significance. Default 0.05.
#' @param epsilon Numeric. Small value added before log. Default 1e-16.
#' @param log_cap Numeric. Max -log10(p) shown. Default is 5.
#' @param point_size Numeric. Size of points. Default is 0.5.
#'
#' @return A ggplot2 object.
#'
#' @details
#' Points are shown as -log10(p), capped at log_cap. If significance columns
#' are provided, point color reflects significance source. Faceted plots show
#' correlation (Pearson r) and a regression line.
#'
#' @export
plot_pvalue_concordance <- function(
  data,
  x_var,
  y_var,
  x_sig_source = NULL,
  y_sig_source = NULL,
  facet_col = NULL,
  facet_levels = NULL,
  x_label = paste("-log10(", x_var, ")"),
  y_label = paste("-log10(", y_var, ")"),
  title = NULL,
  sig_thr = 0.05,
  epsilon = 1e-16,
  log_cap = 5,
  point_size = 0.5
) {

  df <- as.data.frame(data)

  # Error if facet_col is provided but not in the data
  if (!is.null(facet_col) && !(facet_col %in% colnames(df))) {
    stop(sprintf("Column '%s' specified as facet_col not found in data.", facet_col))
  }

  # Optional faceting setup
  if (!is.null(facet_col)) {
    if (is.null(facet_levels)) {
      facet_levels <- head(unique(df[[facet_col]]), 9)
    }
    df <- df[df[[facet_col]] %in% facet_levels, , drop = FALSE]
    df[[facet_col]] <- factor(df[[facet_col]], levels = facet_levels)
  }

  # Compute log-transformed and capped p-values
  df$logp_x <- -log10(df[[x_var]] + epsilon)
  df$logp_y <- -log10(df[[y_var]] + epsilon)
  df$logp_x_capped <- pmin(df$logp_x, log_cap)
  df$logp_y_capped <- pmin(df$logp_y, log_cap)
  df$capped <- factor(ifelse(df$logp_x > log_cap | df$logp_y > log_cap, "Capped", "Uncapped"))

  # Correlation labels and regression lines
  if (!is.null(facet_col)) {
    facet_levels_used <- levels(df[[facet_col]])

    cor_labels <- do.call(rbind, lapply(facet_levels_used, function(facet) {
      sub_df <- df[df[[facet_col]] == facet, ]
      valid_idx <- is.finite(sub_df$logp_x) & is.finite(sub_df$logp_y)
      r_val <- if (sum(valid_idx) >= 2) cor(sub_df$logp_x[valid_idx], sub_df$logp_y[valid_idx]) else NA
      label_text <- paste0("r = ", if (!is.na(r_val)) signif(r_val, 3) else "NA")
      data.frame(label = label_text, facet_value = facet, stringsAsFactors = FALSE)
    }))

    cor_labels[[facet_col]] <- factor(cor_labels$facet_value, levels = levels(df[[facet_col]]))
    cor_labels$facet_value <- NULL

    reg_lines <- do.call(rbind, lapply(facet_levels_used, function(facet) {
      sub_df <- df[df[[facet_col]] == facet, ]
      valid_idx <- is.finite(sub_df$logp_x) & is.finite(sub_df$logp_y)
      fit <- lm(logp_y[valid_idx] ~ logp_x[valid_idx], data = sub_df)
      data.frame(intercept = coef(fit)[1], slope = coef(fit)[2], facet_value = facet)
    }))

    reg_lines[[facet_col]] <- factor(reg_lines$facet_value, levels = levels(df[[facet_col]]))
    reg_lines$facet_value <- NULL

  } else {
    valid_idx <- is.finite(df$logp_x) & is.finite(df$logp_y)
    r_val <- if (sum(valid_idx) >= 2) cor(df$logp_x[valid_idx], df$logp_y[valid_idx]) else NA
    label_text <- paste0("r = ", if (!is.na(r_val)) signif(r_val, 3) else "NA")
    cor_labels <- data.frame(label = label_text)
    fit <- lm(logp_y[valid_idx] ~ logp_x[valid_idx], data = df)
    reg_lines <- data.frame(intercept = coef(fit)[1], slope = coef(fit)[2])
  }

  # Significance coloring
  has_sig <- !is.null(x_sig_source) && !is.null(y_sig_source)
  if (has_sig) {
    df$sig_source <- ifelse(
      df[[x_sig_source]] < sig_thr & df[[y_sig_source]] < sig_thr, "Both",
      ifelse(df[[x_sig_source]] < sig_thr, "Only x",
             ifelse(df[[y_sig_source]] < sig_thr, "Only y", "None"))
    )
  }

  # aes mapping
  aes_mapping <- if (has_sig) {
    aes(
      x = logp_x_capped, 
      y = logp_y_capped, 
      shape = capped, 
      color = sig_source
    )
  } else {
    aes(
      x = logp_x_capped, 
      y = logp_y_capped, 
      shape = capped
    )
  }

  # Plot
  p <- ggplot(
      data = df, 
      mapping = aes_mapping
    ) +
    geom_point(
      size = point_size
    ) +
    geom_abline(
      data = reg_lines, 
      mapping = aes(
        intercept = intercept, 
        slope = slope
      ),
      color = "blue", 
      linetype = "dashed", 
      linewidth = 0.5
    ) +
    scale_shape_manual(
      values = c(
        "Uncapped" = 1, 
        "Capped" = 25
        ), 
      name = "Point status"
    ) 

  if (has_sig) {
    p <- p + scale_color_manual(
      name = "Sig. p-value",
      values = c(
        "None" = "gray80", 
        "Only x" = "#1f77b4", 
        "Only y" = "#ff7f0e", 
        "Both" = "purple"
      ),
      breaks = c(
        "Both", 
        "Only x", 
        "Only y", 
        "None"
      )
    )
  }

  if (!is.null(facet_col)) {
    p <- p + facet_wrap(
        as.formula(
          paste("~", facet_col)
        )
      ) +
      geom_text(
        data = cor_labels, 
        mapping = aes(
          x = -Inf, 
          y = Inf, 
          label = label
        ),
        hjust = -0.1, 
        vjust = 1.1, 
        size = 2, 
        inherit.aes = FALSE
      )
  } else {
    p <- p + geom_text(
        data = cor_labels, 
        mapping = aes(
          x = x, 
          y = y, 
          label = label
        ),
        hjust = -0.1, 
        vjust = 1.1, 
        size = 2, 
        inherit.aes = FALSE
      )
  }

  # Final styling
  p <- p + labs(
      title = title, 
      x = x_label, 
      y = y_label
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend() +
    theme(legend.position = "right")

  return(p)
}

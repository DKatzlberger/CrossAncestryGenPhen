#' Plot correlation differences with optional facets
#'
#' Create a bar plot comparing correlation values across conditions.
#'
#' @param x A data.frame or data.table with required columns.
#' @param cor_cols Two columns with correlation values to compare.
#' @param t_col Column with T statistic values for annotation.
#' @param p_col Column with p-values for annotation.
#' @param facet_col Optional column to use for faceting.
#' @param facet_levels Optional facet values to include in plot.
#' @param x_label Label for the x-axis (default is NULL).
#' @param y_label Label for the y-axis (default is NULL).
#' @param title Optional plot title.
#' @param point_size Point size for annotations.
#'
#' @return A ggplot2 object showing the correlation comparison.
#' 
#' @export
#' @import ggplot2
plot_correlation_difference <- function(
  x,
  cor_cols = c("XR", "YR"),
  t_col,
  p_col,
  facet_col = NULL,
  facet_levels = NULL,
  x_label = NULL,
  y_label = NULL,
  title = NULL,
  point_size = 0.5
) {

  df <- as.data.frame(x)

  # Validate required columns
  required_cols <- c(t_col, p_col, cor_cols)
  if (!is.null(facet_col)) {
    required_cols <- c(required_cols, facet_col)
  }
  if (!all(required_cols %in% names(df))) {
    stop("One or more specified columns do not exist in the data.")
  }

  # Filter by facet levels if specified
  if (!is.null(facet_col) && !is.null(facet_levels)) {
    data <- df[df[[facet_col]] %in% facet_levels, , drop = FALSE]
  }

  # Reshape to long format (base R)
  n <- nrow(df)
  long_data <- data.frame(
    cor_type = rep(cor_cols, each = n),
    cor_value = unlist(df[, cor_cols], use.names = FALSE),
    T_obs = rep(df[[t_col]], times = length(cor_cols)),
    p_value = rep(df[[p_col]], times = length(cor_cols)),
    stringsAsFactors = FALSE
  )
  if (!is.null(facet_col)) {
    long_data[[facet_col]] <- rep(df[[facet_col]], times = length(cor_cols))
  }

  # Construct plot
  p <- ggplot(
    data = long_data, 
    aes(
      x = cor_type, 
      y = cor_value
      )
    ) +
    geom_col(
      width = 0.5
    )

  # Add faceting and annotations if requested
  if (!is.null(facet_col)) {
    p <- p + facet_wrap(
      as.formula(
        paste("~", facet_col)
        )
      )

    # Annotate each facet
    p <- p + geom_text(
      data = long_data[!duplicated(long_data[[facet_col]]), ],
      aes(
        x = -Inf,
        y = Inf,
        label = paste0("T_obs: ", signif(T_obs, 3), "\np-value: ", signif(p_value, 3))
      ),
      hjust = -0.05,
      vjust = 1.2,
      inherit.aes = FALSE,
      size = 2
    )
  } else {
    # Annotate single panel
    p <- p + annotate(
      "text",
      x = 1,
      y = Inf,
      label = paste0("T_obs: ", signif(df[[t_col]][1], 3), "\np-value: ", signif(df[[p_col]][1], 3)),
      hjust = -0.1,
      vjust = 1.2,
      size = 2
    )
  }

  # Final plot decoration
  p <- p +
    labs(
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

#' Plot GDC demographic summaries by a chosen facet
#'
#' Generates demographic summary plots for selected attributes within a
#' clinical or metadata data frame. Each attribute is visualized as either
#' a histogram (for numeric variables) or a bar plot (for categorical variables),
#' faceted by a grouping variable such as \code{POOLED_GENETIC_ANCESTRY}.
#' Optionally, facet labels are only displayed on the last plot to create
#' a more compact layout.
#'
#' @param file_map A data frame containing the metadata or clinical information to be summarized.
#' @param attributes Character vector of column names to plot (required).
#' @param facetting Character string giving the name of the column used for facetting (defaults to \code{"POOLED_GENETIC_ANCESTRY"}).
#' @param facet_last Logical; if \code{TRUE} (default), facet strips are only shown in the last plot.
#'
#' @return A combined \pkg{patchwork} object containing one plot per attribute, arranged side-by-side.
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots
#' @importFrom rlang sym
#' @export
plot_gdc_demographics <- function(
  file_map,
  attributes,
  facetting = "POOLED_GENETIC_ANCESTRY",
  facet_last = TRUE
) {
  ## --- Validation ---
  if (!facetting %in% names(file_map)) stop("data must contain the column specified in facetting.", call. = FALSE)

  ## --- Filter and clean ancestry ---
  df <- as.data.frame(file_map)
  df <- df[!is.na(df[[facetting]]) & df[[facetting]] != "", ]
  df[[facetting]] <- factor(df[[facetting]])

  ## --- Validate attributes ---
  attrs_available <- intersect(attributes, names(df))
  if (length(attrs_available) == 0)stop("None of the specified attributes are present in data", call. = FALSE)

  plots <- list()

  ## --- Generate plots ---
  for (i in seq_along(attrs_available)) {
    col <- attrs_available[i]
    x <- df[[col]]

    if (is.character(x) && all(grepl("^-?[0-9.]+$", na.omit(x)))) {
      x <- as.numeric(x)
      df[[col]] <- x
    }

    df_col <- df[!is.na(x), , drop = FALSE]

    if (is.numeric(x)) {
      p <- ggplot(df_col, aes(x = .data[[col]])) +
        geom_histogram(bins = 30, fill = "grey80", color = "black", linewidth = 0.2) +
        facet_grid(rows = vars(!!sym(facetting)), scales = "free_y") +
        labs(x = col, y = "Frequency") +
        theme_CrossAncestryGenPhen(rotate = 45)
    } else {
      df_col[[col]] <- factor(df_col[[col]])
      p <- ggplot(df_col, aes(x = .data[[col]])) +
        geom_bar(fill = "grey80", color = "black", linewidth = 0.2) +
        facet_grid(rows = vars(!!sym(facetting)), scales = "free_y") +
        labs(x = col, y = "Count") +
        theme_CrossAncestryGenPhen(rotate = 45)
    }

    # Remove facet strips for all but the last plot (if requested)
    if (facet_last && i < length(attrs_available)) {
      p <- p + theme(strip.text.y = element_blank(), strip.background = element_blank())
    }

    plots[[col]] <- p
  }

  ## --- Combine all plots ---
  if (length(plots) == 0) stop("No valid plots could be generated — check your attribute selection.", call. = FALSE)
  combined <- patchwork::wrap_plots(plots, nrow = 1)

  ## --- Return ---
  return(combined)
}

#' Faceted QQ plots for specified genes
#'
#' Plots empirical vs. theoretical quantiles for selected genes, with
#' optional customization of labels and title.
#'
#' @param empirical Matrix of empirical quantiles (quantiles × genes)
#' @param theoretical Matrix of theoretical quantiles (same dimensions)
#' @param method Correlation method: "pearson" or "spearman"
#' @param facet_levels Vector of gene names to plot (must match colnames)
#' @param title Optional plot title
#' @param x_label Optional x-axis label
#' @param y_label Optional y-axis label
#' @param point_size Size of points. Default is 1.
#'
#' @return A ggplot object of QQ plots faceted by gene
#'
#' @import ggplot2
#' @export
plot_qq_correlation <- function(
  empirical,
  theoretical,
  method = c("pearson", "spearman"),
  facet_levels = NULL,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  point_size = 1
){
  # Input checks
  stopifnot(dim(empirical) == dim(theoretical))
  method <- match.arg(method)

  # Vector of gene_names
  feature_names <- colnames(empirical)
  if (is.null(facet_levels)) {
    facet_levels <- feature_names[1:9]
  }
  # Check that features exist
  stopifnot(all(facet_levels %in% feature_names))

  # Subset to selected features
  selected_idx <- match(facet_levels, feature_names)
  selected_emp <- empirical[, selected_idx, drop = FALSE]
  selected_theo <- theoretical[, selected_idx, drop = FALSE]
  selected_names <- feature_names[selected_idx]

  # Correlation for selected features
    # Compute correlation for selected genes
  cors <- sapply(seq_along(selected_idx), function(j) {
    e <- selected_emp[, j]
    t <- selected_theo[, j]
    cor(e, t, method = method)
  })

  df <- data.frame(
    empirical = as.vector(selected_emp),
    theoretical = as.vector(selected_theo),
    feature = factor(
      rep(selected_names, each = nrow(empirical)),
      levels = facet_levels
    )
  )

  method_sub <- paste0(toupper(substr(method, 1, 1)), substr(method, 2, nchar(method)))
  labels <- sprintf("R[%s]==%.2f", method_sub, cors)
  annotation_df <- data.frame(
    feature = factor(facet_levels, levels = facet_levels),
    label = labels
  )

  p <- ggplot(
    data = df,
    mapping = aes(
      x = empirical,
      y = theoretical,
    ) 
  ) +
  geom_point(
    size = point_size
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "gray50",
    linewidth = 0.3
  ) +
  facet_wrap(
    ~feature,
    scales = "free"
  )

  # Annotations
  p  <- p + geom_text(
    data = annotation_df,
    mapping = aes(
      x = -Inf, 
      y = Inf, 
      label = label
    ),
    hjust = -0.1, 
    vjust = 1.1,
    size = 2,
    parse = TRUE,
    inherit.aes = FALSE
  )

  # Final styling
  p <- p + labs(
    title = title,
    x = ifelse(is.null(x_label), "Empirical quantiles", x_label),
    y = ifelse(is.null(y_label), "Theoretical quantiles", y_label)
  ) +
  theme_CrossAncestryGenPhen()

  return(p)
}
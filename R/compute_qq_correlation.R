#' Gene-wise quantile correlations
#'
#' Computes gene-wise correlation between empirical and theoretical quantiles,
#' typically from Q-Q plots, and visualizes the distribution of these correlations
#' as a histogram using ggplot2.
#'
#' @param empirical A numeric matrix of empirical quantiles (quantiles × genes)
#' @param theoretical A numeric matrix of theoretical quantiles (same dimensions as empirical)
#' @param method Character string specifying the correlation method: "pearson" (default) or "spearman"
#' @param title Optional main title for the histogram plot
#' @param x_label Optional label for the x-axis (default depends on method)
#' @param y_label Optional label for the y-axis (default: "Count")
#' @param bins Integer specifying the number of bins for the histogram (default: 50)
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{cors}{A numeric vector of correlation values (length = number of genes)}
#'   \item{method}{The correlation method used}
#'   \item{plot}{A ggplot2 object showing the histogram of correlations}
#' }
#'
#' @import ggplot2
#' @export
compute_qq_correlation <- function(
  empirical,
  theoretical,
  method = c("pearson", "spearman"),
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  bins = 50
){
  ## --- Match method ---
  method <- match.arg(method)

  ## --- Input checks ---
  stopifnot(dim(empirical) == dim(theoretical))


  ## --- Correlations ---
  corrs <- sapply(seq_len(ncol(empirical)), function(j) {
    e <- empirical[, j]
    t <- theoretical[, j]
    cor(e, t, method = method)
  })
  corrs <- setNames(corrs, colnames(empirical))


  ## --- Plot ---
  p <- ggplot(
    data = data.frame(correlation = corrs),
    mapping = aes(
      x = correlation
    )
  ) +
  geom_histogram(
    bins = bins,
    fill = "gray80",
    color = "black",
    linewidth = 0.1
  ) +
  geom_vline(
    aes(
      xintercept = mean(correlation)
    ), 
      color = "red", 
      linewidth = 0.3
  )

  # Final styling
  p <- p + labs(
    title = title,
    x = ifelse(is.null(x_label), paste(paste0(toupper(substr(method, 1, 1)), substr(method, 2, nchar(method))), "correlation"), x_label),
    y = ifelse(is.null(y_label), "Count", y_label),
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_small_legend()

  ## --- Return ---
  return(
    list(
      corrs = corrs,
      plot = p
    )
  )
}
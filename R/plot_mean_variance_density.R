#' Plot Mean-Variance Relationship with Density Overlay
#'
#' This function visualizes the relationship between the mean and variance
#' of log2-transformed expression values across genes. A smooth 2D kernel
#' density estimate is shown as a white-to-blue background, and genes
#' falling in low-density regions are highlighted as black dots.
#'
#' @param X A numeric matrix or data.frame with samples in rows and genes
#'   in columns. Typically the output from limma::voom (i.e., voom(data)$E).
#' @param title Optional title for the plot.
#' @param x_label Optional label for the x-axis.
#' @param y_label Optional label for the y-axis.
#' @param density_grid_size Resolution of the KDE grid (default is 500).
#' @param outlier_quantile Genes below this KDE quantile are shown as outliers (default is 0.75).
#' @param point_size Size of points. Default is 1.
#'
#' @return A ggplot2 object displaying the mean-variance density plot.
#' @import ggplot2
#' @importFrom MASS kde2d
#' @importFrom fields interp.surface
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_mean_variance_density <- function(
  X,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  offset = 0.5,
  density_grid_size = 500, 
  outlier_quantile = 0.75,
  point_size = 1
) {

  # Input validation
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data.frame with samples as rows and genes as columns.")
  }
  

  log2_mean <- apply(X, 2, function(x) log2(mean(x, na.rm = TRUE) + offset))
  log2_variance <- apply(X, 2, function(x) log2(var(x, na.rm = TRUE) + offset))
  df <- data.frame(log2_mean, log2_variance)
  
  # Estimate 2D kernel
  kde <- kde2d(df$log2_mean, df$log2_variance, n = density_grid_size)
  
  # Interpolate density at each gene position
  interp_vals <- interp.surface(obj = kde, loc = df)
  interp_vals[is.na(interp_vals)] <- 0 
  
  # Filter for outliers (e.g., bottom 25% density)
  threshold <- quantile(kde$z, outlier_quantile)
  df_outliers <- df[interp_vals < threshold, ]
  
  # Prepare KDE grid for plotting
  kde_df <- expand.grid(x = kde$x, y = kde$y)
  kde_df$z <- as.vector(kde$z)
  
  # Define custom white-to-blue color scale
  blues <- brewer.pal(9, "Blues")[2:9]
  custom_blues <- c("white", blues)
  
  # Step 7: Generate plot
  p <- ggplot() +
    geom_raster(
      data = kde_df, 
      mapping = aes(
        x = x, 
        y = y, 
        fill = z
      ), 
      interpolate = TRUE
    ) +
    scale_fill_gradientn(
      colors = custom_blues
    ) +
    geom_point(
      data = df_outliers,
      mapping = aes(
        x = log2_mean, 
        y = log2_variance
      ),
      size = point_size
    )

    # Final styling
    p <- p + labs(
      title = title,
      x = ifelse(is.null(x_label), "Log2 Mean expression", x_label),
      y = ifelse(is.null(y_label), " Log2 Variance", y_label)
    ) +
    theme_CrossAncestryGenPhen() +
    theme(legend.position = "none")
  
  return(p)
}

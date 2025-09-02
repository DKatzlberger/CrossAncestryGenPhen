#' t-SNE Cluster Plot
#'
#' Generate a t-SNE plot of combined groups with optional coloring and shaping.
#'
#' @param X Numeric matrix or data frame for the first group.
#' @param Y Numeric matrix or data frame for the second group.
#' @param MX Metadata for X.
#' @param MY Metadata for Y.
#' @param color_var Name of the metadata variable for color (optional).
#' @param shape_var Name of the metadata variable for shape (optional).
#' @param perplexity t-SNE perplexity parameter.
#' @param title Plot title.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param point_size Size of the points.
#' @param seed Optional random seed.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom Rtsne Rtsne
#' @importFrom ggplot2 ggplot aes geom_point labs
plot_tsne_cluster <- function(
  X,
  Y,
  MX,
  MY,
  color_var = NULL,
  shape_var = NULL,
  perplexity = 50,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  point_size = 1,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Combine expression matrices
  X_comb <- rbind(X, Y) 
  M_comb <- rbind(MX, MY)

  # Run t-SNE
  tsne_results <- Rtsne::Rtsne(X_comb, dims = 2, perplexity = perplexity)
  coords <- data.frame(
    TSNE_1 = tsne_results$Y[, 1],
    TSNE_2 = tsne_results$Y[, 2]
  )
  coords <- cbind(coords, M_comb)

  # Build aes() dynamically
  mapping <- aes(x = TSNE_1, y = TSNE_2)
  if (!is.null(color_var)) mapping <- modifyList(mapping, aes(color = .data[[color_var]]))
  if (!is.null(shape_var)) mapping <- modifyList(mapping, aes(shape = .data[[shape_var]]))

  # Plot
  p <- ggplot(
    data = coords, 
    mapping = mapping
    ) +
    geom_point(
      size = point_size
    ) +
    labs(
      title = title,
      x = ifelse(is.null(x_label), "TSNE1", x_label),
      y = ifelse(is.null(y_label), "TSNE2", y_label),
      color = color_var,
      shape = shape_var
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  return(p)
}

#' PCA Cluster Plot
#'
#' Generate a PCA plot of combined groups with optional coloring and shaping.
#' Axis labels always include variance explained for the selected PCs.
#'
#' @param X Numeric matrix or data frame for the first group.
#' @param Y Numeric matrix or data frame for the second group.
#' @param MX Metadata for X.
#' @param MY Metadata for Y.
#' @param color_var Name of the metadata variable for color (optional).
#' @param shape_var Name of the metadata variable for shape (optional).
#' @param pcs Integer vector of length 2 indicating which PCs to plot (default c(1,2)).
#' @param title Plot title.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param point_size Size of the points.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point labs
plot_pca_cluster <- function(
  X,
  Y,
  MX,
  MY,
  color_var = NULL,
  shape_var = NULL,
  pcs = c(1, 2),
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  point_size = 1
) {
  # Combine expression matrices & metadata
  X_comb <- rbind(X, Y)
  M_comb <- rbind(MX, MY)

  # PCA (center & scale are standard for clustering/visualization)
  pca <- stats::prcomp(X_comb, center = TRUE, scale. = TRUE)

  # Variance explained
  var_exp <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)

  # Validate PCs
  stopifnot(length(pcs) == 2, all(pcs >= 1), max(pcs) <= ncol(pca$x))
  x_pc <- pcs[1]
  y_pc <- pcs[2]

  # Prepare coordinates + metadata
  coords <- data.frame(
    PC_X = pca$x[, x_pc],
    PC_Y = pca$x[, y_pc],
    M_comb,
    check.names = FALSE
  )

  # Build aes() dynamically
  mapping <- ggplot2::aes(x = .data[["PC_X"]], y = .data[["PC_Y"]])
  if (!is.null(color_var)) mapping <- modifyList(mapping, ggplot2::aes(color = .data[[color_var]]))
  if (!is.null(shape_var)) mapping <- modifyList(mapping, ggplot2::aes(shape = .data[[shape_var]]))

  # Axis labels with variance explained
  xl <- sprintf("PC%d (%.1f%%)", x_pc, 100 * var_exp[x_pc])
  yl <- sprintf("PC%d (%.1f%%)", y_pc, 100 * var_exp[y_pc])

  # Plot
  p <- ggplot2::ggplot(
      data = coords,
      mapping = mapping
    ) +
    ggplot2::geom_point(
      size = point_size
    ) +
    ggplot2::labs(
      title = title,
      x = ifelse(is.null(x_label), xl, x_label),
      y = ifelse(is.null(y_label), yl, y_label),
      color = color_var,
      shape = shape_var
    ) +
    theme_CrossAncestryGenPhen()

  return(p)
}

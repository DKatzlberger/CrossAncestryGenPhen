#' Expression Heatmap with Ancestry and Group Split
#'
#' Generates a z-scored expression heatmap with samples split by combined ancestry and group,
#' and annotations for ancestry and group using ComplexHeatmap.
#'
#' @param X Expression matrix for the first ancestry.
#' @param Y Expression matrix for the second ancestry.
#' @param MX Metadata for X.
#' @param MY Metadata for Y.
#' @param a_col Name of the column in metadata that defines ancestry.
#' @param g_col Name of the column in metadata that defines groups.
#' @param features Optional vector of feature names to plot. Defaults to the first 9 features.
#' @param title Optional heatmap title.
#' @param x_label Optional label for the color legend.
#' @param y_label Optional label for the row axis.
#' @param file Optional file path to save the plot (PDF, PNG, JPEG, or TIFF).
#' @param width Plot width in inches (default = 8).
#' @param height Plot height in inches (default = 8).
#'
#' @return Invisibly returns NULL. Called for its side effect of drawing or saving a heatmap.
#' @export
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom circlize colorRamp2
plot_expression_heatmap <- function(
  X,
  Y,
  MX,
  MY,
  a_col,
  g_col,
  features = NULL,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  file = NULL,
  width = 8,
  height = 8
){

  X <- rbind(X, Y)
  M <- rbind(MX, MY)

  # Select features
  if (is.null(features)) {
    features <- colnames(X)
    features <- sort(features)
    features <- features[1:min(9, ncol(X))]
  } 
  features <- intersect(features, colnames(X))
  if (length(features) == 0) stop("None of the selected features are found.")

  # Z-score transformation
  mat <- t(scale(X[, features]))

  z_min <- min(mat, na.rm = TRUE)
  z_max <- max(mat, na.rm = TRUE)
  z_lim <- max(abs(c(z_min, z_max)))


  # Make sure a_col and g_col are in M
  if (!(a_col %in% colnames(M))) stop("a_col not found in metadata.")
  if (!(g_col %in% colnames(M))) stop("g_col not found in metadata.")

  a_factor <- factor(M[[a_col]])
  a_levels <- levels(a_factor)
  a_1 <- levels(a_factor)[1]
  a_2 <- levels(a_factor)[2]

  g_factor <- factor(M[[g_col]])
  g_levels <- levels(g_factor)
  g_1 <- levels(g_factor)[1]
  g_2 <- levels(g_factor)[2]

  # Colors
  a_colors <- c("#3399ff", "#ffe76d")  # Green, Purple
  g_colors <- c("#6a3d9a", "#1f9e89")  # Brown-Gold, Teal

  # Named list
  col_list <- list(
    setNames(a_colors, a_levels),
    setNames(g_colors, g_levels)
  )
  names(col_list) <- c(a_col, g_col)

  # Combine factors into named list
  anno_list <- setNames(list(a_factor, g_factor), c(a_col, g_col))

  # Create annotation
  top_anno <- do.call(
    ComplexHeatmap::HeatmapAnnotation, 
    c(
      anno_list,
      list(
        col = col_list,
        border = TRUE,
        show_annotation_name = FALSE,
        annotation_legend_param = list(direction = "horizontal")
      )
    )
  )

  # Hierarchical split: ancestry then condition
  column_split <- data.frame(
    Ancestry = a_factor,
    Condition = g_factor
  )

  # Heatmap
  ht <- Heatmap(
    mat,
    col = colorRamp2(
      c(-z_lim, 0, z_lim), 
      c("blue", "white", "red")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    column_title = title,
    row_title = ifelse(is.null(y_label), "feature", y_label),
    name = ifelse(is.null(x_label), "z-score", x_label),
    row_title_side = "left",
    row_names_side = "left",
    column_title_side = "top",
    column_split = column_split,
    top_annotation = top_anno,
    heatmap_legend_param = list(direction = "horizontal"),
  )

  if (!is.null(file)) {
    ext <- tools::file_ext(file)
    if (ext == "pdf") {
      pdf(file, width = width, height = height)
      draw(
        ht,
        heatmap_legend_side = "bottom",
        annotation_legend_side = "bottom",
        merge_legend = TRUE
        )
      dev.off()
    } else if (ext %in% c("png", "jpeg", "tiff")) {
      do.call(ext, list(filename = file, width = width, height = height, units = "in", res = 300))
      draw(
        ht,  
        heatmap_legend_side = "bottom",
        annotation_legend_side = "bottom",
        merge_legend = TRUE
      )
      dev.off()
    } else {
      stop("Unsupported file extension. Use pdf, png, jpeg, or tiff.")
    }
  } else {
    draw(
      ht,  
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      merge_legend = TRUE
    )
  }

  invisible(NULL)
}
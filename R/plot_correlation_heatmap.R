#' ComplexHeatmap of Correlation Matrix Across Iterations
#'
#' Computes correlation between values (e.g., p-values, t-statistics) across iterations per feature,
#' reshapes the data into wide format, calculates a correlation matrix, and visualizes it as a clustered heatmap.
#'
#' @param x A data frame with columns: `feature`, `iteration`, and a numeric value column (e.g., `p_value`, `T_obs`).
#' @param value_col A string specifying the column name in `x` to use for correlation (e.g., `"p_value"` or `"T_obs"`).
#' @param method A string specifying the correlation method to use. One of `"pearson"` (default), `"spearman"`, or `"kendall"`.
#' @param row_names Optional string to display as the row title in the heatmap.
#' @param title Optional string for the main heatmap title.
#' @param file Optional file path to save the plot (e.g., `"plot.pdf"` or `"plot.png"`). If `NULL`, the plot is shown interactively.
#' @param width Width of the plot in inches (default = 8).
#' @param height Height of the plot in inches (default = 8).
#'
#' @return Invisibly returns `NULL`. The function is called for its side effect of drawing or saving a heatmap.
#'
#' @details
#' This function is useful for visualizing the stability or similarity of a given statistic (like p-values or test statistics)
#' across multiple iterations or resampling runs. It internally uses `compute_correlation_matrix()` to compute the correlation matrix.
#'
#' The upper triangle of the matrix is used to calculate a summary mean correlation, which is displayed above the heatmap.
#'
#' @import data.table
#' @import ComplexHeatmap
#' @import circlize
#' @import grid
#' @export

plot_correlation_heatmap <- function(
  x,
  value_col,
  method = "pearson",
  row_names = NULL,
  title = NULL,
  file = NULL,
  width = 8,
  height = 8
) {

  mat <- compute_correlation_matrix(x, value_col, method = method)
  diag(mat) <- NA 

  # Mean correlation across iterations
  mask <- upper.tri(mat, diag = FALSE)
  mean <- mean(mat[mask], na.rm = TRUE)

  # Annotation for heatmap
  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    anno = function(index) {
      grid::grid.text(
        label = paste0("Mean Pearson coefficient: ", round(mean, 3)),
        just = "left",
        x = unit(0, "npc"),
        y = unit(0.5, "npc"),
        gp = grid::gpar(fontface = "italic")
      )
    },
    show_annotation_name = FALSE,
    which = "column",
    height = unit(1.5, "lines")
  )
  

  ht <- Heatmap(
    mat,
    name = "Pearson",
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    row_title = if (!is.null(row_names)) as.character(row_names) else NULL,
    column_title = if (!is.null(title)) as.character(title) else NULL,
    row_title_side = "right",
    column_title_side = "top",
    column_names_rot = 0,
    top_annotation = top_anno,
    heatmap_legend_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
  )

  if (!is.null(file)) {
    ext <- tools::file_ext(file)
    if (ext == "pdf") {
      pdf(file, width = width, height = height)
      draw(ht, heatmap_legend_side = "right")
      dev.off()
    } else if (ext %in% c("png", "jpeg", "tiff")) {
      do.call(ext, list(filename = file, width = width, height = height, units = "in", res = 300))
      draw(ht, heatmap_legend_side = "right")
      dev.off()
    } else {
      stop("Unsupported file extension. Use pdf, png, jpeg, or tiff.")
    }
  } else {
    draw(ht, heatmap_legend_side = "right")
  }

  invisible(NULL)
}

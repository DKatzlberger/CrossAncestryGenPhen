#' ComplexHeatmap of P-value Correlations Across Iterations
#'
#' Computes correlation between p-values per iteration and plots with clustering.
#'
#' @param x Data frame with `feature`, `iteration`, and `p_value` columns.
#' @param file Optional path to save the plot (e.g., "plot.pdf" or "plot.png"). If NULL, shows interactively.
#' @param width Plot width in inches (default = 8).
#' @param height Plot height in inches (default = 8).
#'
#' @return Draws the heatmap; optionally saves it.
#' @import data.table
#' @import ComplexHeatmap
#' @import circlize
#' @export
plot_repeated_pvalue_correlation <- function(
  x,
  file = NULL,
  width = 8,
  height = 8
) {

  dt <- as.data.table(x)
  wide <- dcast(dt, feature ~ iteration, value.var = "p_value")
  pval_mat <- as.matrix(wide[, -1, with = FALSE])
  cor_mat <- cor(pval_mat, use = "pairwise.complete.obs", method = "pearson")
  diag(cor_mat) <- NA

  ht <- Heatmap(
    cor_mat,
    name = "Pearson",
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    show_row_names = FALSE,        
    show_column_names = FALSE,
    column_title = "Correlation of P-values Across Iterations",
    heatmap_legend_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
  )

  if (!is.null(file)) {
    ext <- tools::file_ext(file)
    if (ext == "pdf") {
      pdf(file, width = width, height = height)
      draw(ht)
      dev.off()
    } else if (ext %in% c("png", "jpeg", "tiff")) {
      do.call(ext, list(filename = file, width = width, height = height, units = "in", res = 300))
      draw(ht)
      dev.off()
    } else {
      stop("Unsupported file extension. Use pdf, png, jpeg, or tiff.")
    }
  } else {
    draw(ht)
  }
}

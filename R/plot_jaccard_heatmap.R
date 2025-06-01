#' Plot Jaccard Heatmap of Sample Reuse Across Iterations
#'
#' Computes and visualizes the Jaccard similarity of sample usage across iterations
#' for a given role (e.g., "test", "train", or "inference").
#'
#' @param id_usage A data.frame with columns \code{sample_id}, \code{role}, and \code{iteration}.
#' @param role A string, one of \code{"test"}, \code{"train"}, or \code{"inference"}.
#' @param title Optional title for the heatmap.
#' @param file Optional filename to save the heatmap (pdf, png, jpeg, or tiff).
#' @param width, height Width and height of the saved file (in inches).
#'
#' @return Invisibly returns the Jaccard matrix (numeric matrix).
#'
#' @export
plot_jaccard_heatmap <- function(
  id_usage,
  role = "test",
  row_names = NULL,
  title = NULL,
  file = NULL,
  width = 8,
  height = 8
) {

  mat <- compute_jaccard_matrix(id_usage, role = role)
  diag(mat) <- NA 
  
  # Mean jaccard similarity across iterations
  mask <- upper.tri(mat, diag = FALSE)
  mean <- mean(mat[mask], na.rm = TRUE)

  # Annotation for heatmap
  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    anno = function(index) {
      grid::grid.text(
        label = paste0("Mean Jaccard index: ", round(mean, 3)),
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

  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "Jaccard",
    col = circlize::colorRamp2(c(0, 1), c("white", "darkgreen")),
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
    heatmap_legend_param = list(
      at = c(0, 0.5, 1),
      labels = c("0", "0.5", "1")
    )
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

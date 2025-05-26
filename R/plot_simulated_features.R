#' Plot Simulated Features by Ancestry and Group (base R + data.table)
#'
#' Visualize expression of selected genes from simulated data, grouped by ancestry and filled by group.
#'
#' @param sim_data Output from `simulate_DEG_signal()` or `simulate_NB_counts()`.
#' @param features Optional character vector of gene names to plot. If NULL, samples `n` genes.
#' @param a_col Column in `sample_info` used for x-axis grouping (e.g., "ancestry").
#' @param g_col Column in `sample_info` used for fill color (e.g., "group").
#' @param zscore Logical. If TRUE, z-score normalize each gene across samples.
#' @param n Number of genes to sample if `features` is NULL.
#' @param seed Optional random seed.
#'
#' @return A ggplot2 object showing boxplots for selected features.
#' @export
plot_simulated_features <- function(
  sim_data,
  features = NULL,
  a_col = "ancestry",
  g_col = "group",
  zscore = TRUE,
  n = 12,
  seed = NULL
) {
  library(data.table)
  library(ggplot2)

  if (!is.null(seed)) set.seed(seed)

  counts <- sim_data$counts
  info <- sim_data$sample_info

  all_genes <- rownames(counts)
  if (is.null(features)) {
    features <- sample(all_genes, n)
  }

  # Subset and optionally z-score
  dat <- counts[features, , drop = FALSE]
  if (zscore) {
    dat <- t(scale(t(dat)))
  }

  # Create long-format data.table
  expr_dt <- as.data.table(dat)
  expr_dt[, feature := rownames(dat)]
  expr_long <- melt(expr_dt, id.vars = "feature", variable.name = "sample", value.name = "expression")

  # Add idx to both tables for reliable join
  info$idx <- seq_len(nrow(info))
  expr_long[, idx := match(sample, colnames(dat))]

  # Merge with sample info
  info_dt <- as.data.table(info)
  merged_dt <- merge(expr_long, info_dt, by = "idx", all.x = TRUE)

  # Plot
  ggplot(merged_dt, aes_string(x = a_col, y = "expression", fill = g_col)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
    facet_wrap(~ feature, scales = "free_y", ncol = 4) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      x = a_col,
      y = if (zscore) "Z-scored expression" else "Expression",
      fill = g_col,
      title = "Simulated Gene Expression by Ancestry and Group"
    )
}

#' Simulate two populations with synthetic bulk RNA-seq count data (compcodeR)
#'
#' Uses compcodeR to simulate two populations with balanced group sizes and
#' independent differential expression. Intended log2 fold-change is 
#' recorded for each population.
#'
#' @param n_vars Total number of features (genes).
#' @param n_degs_X Number of DEGs in population X.
#' @param n_degs_Y Number of DEGs in population Y.
#' @param log2FC_X Intended log2 fold-change for pop X.
#' @param log2FC_Y Intended log2 fold-change for pop Y.
#' @param n_X Number of samples per condition for pop X.
#' @param n_Y Number of samples per condition for pop Y.
#' @param seed Optional seed for reproducibility.
#'
#' @return A list with:
#' \describe{
#'   \item{X}{Count matrix for population X}
#'   \item{Y}{Count matrix for population Y}
#'   \item{MX}{Sample metadata for X}
#'   \item{MY}{Sample metadata for Y}
#'   \item{fX}{Feature truth table for X}
#'   \item{fY}{Feature truth table for Y}
#'   \item{summary}{Population and condition summary tables}
#' }
#' @import compcodeR
#' @export
sim_imbalanced_groups <- function(
  n_vars,
  n_degs_X,
  n_degs_Y,
  log2FC_X,
  log2FC_Y,
  n_X,
  n_Y,
  seed = NULL
) {

  # Helper to generate synthetic data
  sim_pop <- function(
    a_name,
    n_samples,
    n_vars,
    n_degs,
    log2FC,
    seed = NULL
  ) {

    if (!is.null(seed)) set.seed(seed)

    # Use generateSyntheticData from compcodeR
    sim <- compcodeR::generateSyntheticData(
      dataset = a_name,
      n.vars = n_vars,
      samples.per.cond = n_samples,
      n.diffexp = n_degs,
      seqdepth = 1e7,
      effect.size = 2^log2FC,
      fraction.upregulated = 0.5,
      filter.threshold.total = 0,
      filter.threshold.mediancpm = 0
    )

    # Counts
    counts <- sim@count.matrix
    colnames(counts) <- paste0(a_name, "_", seq_len(ncol(counts)))

    # Meta
    conditions <- sim@sample.annotations$condition
    conditions <- factor(conditions, levels = c(1, 2), labels = c("H", "D"))

    meta <- data.frame(
      id = colnames(counts),
      condition = conditions,
      ancestry = a_name,
      stringsAsFactors = FALSE
    )

    # Features
    is_H <- meta$condition == "H"
    is_D <- meta$condition == "D"
    mean_H <- rowMeans(counts[, is_H, drop = FALSE] + 1)
    mean_D <- rowMeans(counts[, is_D, drop = FALSE] + 1)
    observed_log2FC <- log2(mean_D / mean_H)

    features <- data.frame(
      feature = rownames(counts),
      is_DE = sim@variable.annotations$differential.expression,
      intended_log2FC = log2FC,
      true_log2FC = sim@variable.annotations$truelog2foldchanges,
      observed_log2FC = observed_log2FC,
      ancestry = a_name
    )

    list(
      counts = counts,
      meta = meta,
      features = features
    )
  }

  # Simulation step
  X <- sim_pop("X", n_X, n_vars, n_degs_X, log2FC_X, if (!is.null(seed)) seed else NULL)
  Y <- sim_pop("Y", n_Y, n_vars, n_degs_Y, log2FC_Y, if (!is.null(seed)) seed + 1000 else NULL)

  fX <- X$features
  fY <- Y$features

  # Interaction truth
  interaction <- data.frame(
    feature = fX$feature,
    is_DE = as.numeric(xor(fX$is_DE, fY$is_DE)),
    intended_log2FC = log2FC_Y - log2FC_X,
    true_log2FC = fY$true_log2FC - fX$true_log2FC,
    observed_log2FC = fY$observed_log2FC - fX$observed_log2FC,
    ancestry = "Interaction"
  )

  # Summary statistic
  # Demographics (Condition/Population)
  table_X <- table(X$meta$condition)
  table_Y <- table(Y$meta$condition)

  ancestry_summary <- data.frame(
    ancestry = c("X", "Y"),
    total_samples = c(sum(table_X), sum(table_Y)),
    ratio = round(
      c(
        (sum(table_X) / sum(table_Y)), 
        1
      ),
      2
    )
  )

  condition_summary <- data.frame(
    ancestry = c("X", "Y"),
    H = c(table_X["H"], table_Y["H"]),
    D = c(table_X["D"], table_Y["D"]),
    ratio = round(
      c(
        (table_X["H"] / table_X["D"]), 
        (table_Y["H"] / table_Y["D"])
      ),
      2
    )
  )

  # Feature level
  feature_summary <- rbind(
    data.frame(
      ancestry = "X",
      total_features = nrow(fX),
      total_DEGs = sum(fX$is_DE == 1, na.rm = TRUE),
      total_non_DEGs = sum(fX$is_DE == 0, na.rm = TRUE),
      pct_DEGs = mean(fX$is_DE == 1, na.rm = TRUE) * 100
    ),
    data.frame(
      ancestry = "Y",
      total_features = nrow(fY),
      total_DEGs = sum(fY$is_DE == 1, na.rm = TRUE),
      total_non_DEGs = sum(fY$is_DE == 0, na.rm = TRUE),
      pct_DEGs = mean(fY$is_DE == 1, na.rm = TRUE) * 100
    )
  )

  # Interaction summary
  interaction_summary <- data.frame(
    ancestry = "Interaction",
    total_features = nrow(interaction),
    total_DEGs = sum(interaction$is_DE == 1, na.rm = TRUE),
    total_non_DEGs = sum(interaction$is_DE == 0, na.rm = TRUE),
    pct_DEGs = mean(interaction$is_DE) * 100
  )


  # Bundle all
  summary <- list(
    ancestry = ancestry_summary,
    condition = condition_summary,
    feature = feature_summary,
    interaction = interaction_summary
  )

  return(
    list(
      X = t(X$counts),
      Y = t(Y$counts),
      MX = X$meta,
      MY = Y$meta,
      fX = X$features,
      fY = Y$features,
      interaction = interaction,
      summary = summary
    )
  )
}

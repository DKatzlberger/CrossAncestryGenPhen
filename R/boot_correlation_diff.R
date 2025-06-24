#' Bootstrap-Based Correlation Difference Test (Unstratified)
#'
#' Estimates uncertainty in the correlation difference between ancestries X and Y
#' (relative to reference R), using non-stratified bootstrap resampling within ancestries.
#'
#' @param X Expression matrix for ancestry X (samples x genes)
#' @param Y Expression matrix for ancestry Y
#' @param R Expression matrix for reference ancestry
#' @param MX Metadata for X
#' @param MY Metadata for Y
#' @param MR Metadata for R
#' @param g_col Column in metadata specifying condition/group (factor with 2 levels)
#' @param method Correlation method: "pearson" (default) or "spearman"
#' @param B Number of bootstrap samples
#' @param seed Optional seed for reproducibility
#' @param alpha Significance level (e.g., 0.05 for 95\% CI)
#'
#' @export
boot_correlation_diff <- function(
  X,
  Y,
  R,
  MX,
  MY,
  MR,
  g_col,
  method = c("pearson", "spearman"),
  B = 1000,
  seed = NULL,
  alpha = 0.05
) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)

  # Extract group labels
  g_X <- MX[[g_col]]
  g_Y <- MY[[g_col]]
  g_R <- MR[[g_col]]
  g1 <- levels(factor(g_X))[1]
  g2 <- levels(factor(g_X))[2]

  # Fixed reference contrast
  delta_R <- colMeans(R[g_R == g1, , drop = FALSE]) - colMeans(R[g_R == g2, , drop = FALSE])

  # Observed deltas and correlation
  delta_X <- colMeans(X[g_X == g1, , drop = FALSE]) - colMeans(X[g_X == g2, , drop = FALSE])
  delta_Y <- colMeans(Y[g_Y == g1, , drop = FALSE]) - colMeans(Y[g_Y == g2, , drop = FALSE])
  cor_XR <- cor(delta_X, delta_R, use = "complete.obs", method = method)
  cor_YR <- cor(delta_Y, delta_R, use = "complete.obs", method = method)
  T_obs <- cor_YR - cor_XR

  # Bootstrap resampling without preserving group structure
  T_boot <- rep(NA_real_, B)
  B_used <- 0

  for (b in seq_len(B)) {
    # Sample rows from X and Y (non-stratified)
    Xb_idx <- sample(nrow(X), replace = TRUE)  
    Yb_idx <- sample(nrow(Y), replace = TRUE)

    Xb <- X[Xb_idx, , drop = FALSE]
    Yb <- Y[Yb_idx, , drop = FALSE]

    g_Xb <- g_X[Xb_idx]
    g_Yb <- g_Y[Yb_idx]

    # Check if both groups are present in both samples
    if (!all(c(g1, g2) %in% g_Xb) || !all(c(g1, g2) %in% g_Yb)) next

    # Compute resampled condition effects
    dX <- colMeans(Xb[g_Xb == g1, , drop = FALSE]) - colMeans(Xb[g_Xb == g2, , drop = FALSE])
    dY <- colMeans(Yb[g_Yb == g1, , drop = FALSE]) - colMeans(Yb[g_Yb == g2, , drop = FALSE])

    # Compute correlation difference
    c1 <- cor(dX, delta_R, use = "complete.obs", method = method)
    c2 <- cor(dY, delta_R, use = "complete.obs", method = method)
    if (is.na(c1) || is.na(c2)) next

    T_boot[B_used + 1] <- c2 - c1
    B_used <- B_used + 1
  }

  T_boot <- T_boot[seq_len(B_used)]
  if (B_used == 0) stop("All bootstrap samples failed.")

  # Summary statistics
  SE <- sd(T_boot, na.rm = TRUE)
  CI <- quantile(T_boot, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

  summary_stats <- data.frame(
    feature = "Global",
    T_obs = T_obs,
    cor_XR = cor_XR,
    cor_YR = cor_YR,
    SE = SE,
    CI_lower = CI[1],
    CI_upper = CI[2],
    p_value = NA_real_,
    row.names = NULL
  )

  return(list(
    summary_stats = summary_stats,
    T_replicates = matrix(T_boot, ncol = 1, dimnames = list(NULL, "diff_cor")),
    B_used_replicates = B_used
  ))
}

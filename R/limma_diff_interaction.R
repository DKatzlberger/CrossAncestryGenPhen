#' Internal: Limma-Based Interaction Test
#'
#' Computes gene-wise interaction statistics using limma for the difference in group effects
#' between ancestries (ancestry X minus ancestry Y), matching the logic of
#' \code{perm_diff_interaction}.
#'
#' @param X Expression matrix for ancestry X. Rows = samples, columns = genes.
#' @param Y Expression matrix for ancestry Y. Rows = samples, columns = genes.
#' @param MX Metadata for X. Must include group and ancestry columns.
#' @param MY Metadata for Y. Must include group and ancestry columns.
#' @param g_col Name of the column indicating group.
#' @param a_col Name of the column indicating ancestry.
#'
#' @return A list with:
#' \describe{
#'   \item{summary_stats}{A data.frame with gene-level interaction coefficients, p-values, FDR, and confidence intervals.}
#'   \item{fit}{The full limma model fit.}
#'   \item{group_levels}{Factor levels used for group variable.}
#'   \item{coef_tested}{Name of the interaction coefficient extracted.}
#' }
#' @keywords internal
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats qt
#' @noRd
limma_diff_interaction <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col
) {
  stopifnot(is.matrix(X), is.matrix(Y))
  stopifnot(ncol(X) == ncol(Y))

  # Extract and validate group/ancestry
  g_X <- MX[[g_col]]
  g_Y <- MY[[g_col]]
  validate_groups(g_X, g_Y)
  g1 <- levels(g_X)[1]
  g2 <- levels(g_X)[2]

  a_X <- unique(MX[[a_col]])
  a_Y <- unique(MY[[a_col]])
  validate_ancestry(a_X, a_Y)

  # Combine expression and metadata
  XY <- rbind(X, Y)
  MXY <- rbind(MX, MY)

  # Factor setup
  group <- factor(MXY[[g_col]], levels = c(g1, g2)) 
  ancestry <- factor(MXY[[a_col]], levels = c(a_X, a_Y))  

  # Design matrix with interaction
  design <- model.matrix(~ group * ancestry)

  # Fit linear model
  fit <- limma::lmFit(t(XY), design)
  fit <- limma::eBayes(fit)

  # Identify interaction coefficient
  coef_name <- grep("^group.*:ancestry", colnames(design), value = TRUE)
  if (length(coef_name) != 1) {
    stop("Interaction term not uniquely identified.")
  }

  # Extract statistics
  res <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  # Standard error
  res$SE <- res$logFC / res$t

  summary_stats <- data.frame(
    feature = rownames(res),
    T_obs = res$logFC,
    SE = res$SE,
    p_value = res$P.Value,
    p_adj = res$adj.P.Val,
    row.names = NULL
  )

  return(list(
    summary_stats = summary_stats,
    fit = fit,
    group_levels = levels(group),
    ancestry_levels = levels(ancestry),
    coef = coef_name
    )
  )
}

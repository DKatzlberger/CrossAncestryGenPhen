#' Summarize subsets results by feature/coefficient with multiple aggregation methods
#'
#' @param stats A data frame with columns: `feature`, `T_obs`, `p_value`, `p_adj`, `ave_expr`.
#' @param method Aggregation method to highlight/return separately. Options: "mean", "cct", "fisher", "bonferroni".
#' @param by Additional grouping variables besides `feature` (character vector).
#'
#' @importFrom data.table as.data.table
#' @importFrom data.table .SD
#' @importFrom stats p.adjust
#'
#' @export
summarize_subsets <- function(
  stats,
  method = c("mean", "cct", "fisher", "bonferroni"),
  by = NULL
) {
  method <- match.arg(method)


  ## --- Check input ---
  required <- c(
    "coef_id", "coef_type", "contrast", "g_1", "g_2", "a_1", "a_2", 
    "feature", "T_obs", "SE", "p_value", "p_adj", "ave_expr"
  )
  
  missing  <- setdiff(required, colnames(stats))
  if (length(missing) > 0) {
    stop("[summarize_subsets] Missing required column(s): ", paste(missing, collapse = " "))
  }


  ## --- Check grouping vars ---
  if (!is.null(by) && !all(by %in% colnames(stats))) {
    stop("[summarize_subsets] Grouping column(s) not found: ",
         paste(setdiff(by, colnames(stats)), collapse = " "))
  }
  per_feature_grouping <- unique(c("feature", by))


  ## --- Helper functions ---
  mean_pval <- function(pvec) {
    pvec <- pvec[!is.na(pvec)]
    if (!length(pvec)) return(NA_real_)
    mean(pvec)
  }

  cct_pval <- function(pvec) {
    pvec <- pvec[!is.na(pvec)]
    if (!length(pvec)) return(NA_real_)
    pvec <- pmin(pmax(pvec, 1e-15), 1 - 1e-15)
    stat <- mean(tan((0.5 - pvec) * pi))
    p <- 0.5 - atan(stat) / pi
    max(min(p, 1), 0)
  }

  fisher_pval <- function(pvec) {
    pvec <- pvec[!is.na(pvec)]
    if (!length(pvec)) return(NA_real_)
    stat <- -2 * sum(log(pvec))
    df   <- 2 * length(pvec)
    pchisq(stat, df = df, lower.tail = FALSE)
  }

  bonferroni_pval <- function(pvec) {
    pvec <- pvec[!is.na(pvec)]
    if (!length(pvec)) return(NA_real_)
    min(1, min(pvec) * length(pvec))
  }

  aggregators <- list(
    mean       = mean_pval,
    cct        = cct_pval,
    fisher     = fisher_pval,
    bonferroni = bonferroni_pval
  )


  ## --- Data.table prep ---
  dt <- as.data.table(stats)
  keep_cols <- setdiff(colnames(dt), c("iteration","T_obs", "SE", "p_value","p_adj","ave_expr", per_feature_grouping))


  ## --- Run each aggregation method ---
  all_results <- lapply(
    names(aggregators), function(m) {
      agg_dt <- dt[, {
        T_obs    <- mean(T_obs, na.rm = TRUE)
        SE       <- mean (SE, na.rm = TRUE)
        ave_expr <- mean(ave_expr, na.rm = TRUE)

        res <- list(
          T_obs    = T_obs,
          SE       = SE,
          ave_expr = ave_expr,
          p_value  = aggregators[[m]](p_value)
        )

        for (col in keep_cols) res[[col]] <- .SD[[col]][1]
        res
      }, by = per_feature_grouping]

      ## Adjust p-values across features per contrast (loop over agg methods)
      per_contrast_grouping <- setdiff(per_feature_grouping, "feature")
      agg_dt[, p_adj := p.adjust(p_value, method = "BH"), by = per_contrast_grouping]

      ## Explicit forced order (same as limma_interaction_effect)
      out_cols <- c(
        "coef_id", "coef_type", "contrast", "g_1", "g_2", "a_1", "a_2",
        "feature", "T_obs", "SE", "p_value", "p_adj", "ave_expr"
      )

      ## Keep only those columns that actually exist
      out_cols <- out_cols[out_cols %in% colnames(agg_dt)]

      as.data.frame(agg_dt[, ..out_cols])
    }
  )

  ## --- Add method name ---
  names(all_results) <- names(aggregators)


  ## --- Return both all + selected ---
  out <- list(
    all_method = all_results,
    sel_method = all_results[[method]]
  )

  
  ## --- Return ---
  return(out)
}

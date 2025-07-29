#' Summarize subsets results by feature, with Cauchy-combined p-value
#'
#' @param data A data frame with columns: `feature`, `T_obs`, `p_value`, `p_adj`, `ave_expr`, `iteration`.
#' @param iter_col Name of the iteration column (string).
#'
#' @export
summarize_subsets <- function(
  data,
  iter_col
) { 

  # Harmonic mean helper
  hm_rank <- function(
    ranks, 
    na.rm = FALSE
  ) {
    if (na.rm) ranks <- ranks[!is.na(ranks)]
    ranks <- ranks[ranks > 0]
    if (length(ranks) == 0) return(NA_real_)
    length(ranks) / sum(1 / ranks)
  }

  # Cauchy combination helper
  cct_pval <- function(
    pvec,
    na.rm = FALSE,
    epsilon = 1e-15
  ) {
    # Remove NA pvals and make them finite
    if (na.rm) (pvec <- pvec[!is.na(pvec)])
    if (length(pvec) == 0) return(NA_real_)
    pvec <- pmin(pmax(pvec, epsilon), 1 - epsilon)

    # Cauchy statistic
    tstat <- mean(tan((0.5 - pvec) * pi))
    p_comb <- 1 - pcauchy(tstat)
    p_comb <- min(max(p_comb, 0), 1) 

    return(p_comb)
  }

  # Add ranks
  data$rank <- ave(
    data$p_value,
    data[[iter_col]],
    FUN = function(x) rank(x, ties.method = "average")
  )
  
  # Group by feature
  agg_list <- by(data, data$feature, function(sub) {

      T_obs <- mean(sub$T_obs, na.rm = TRUE)
      p_value <- mean(sub$p_value, na.rm = TRUE)
      cct_value <- cct_pval(sub$p_value, na.rm = TRUE)
      prob_sig <- mean(sub$p_adj < 0.05, na.rm = TRUE)
      ave_expr <- mean(sub$ave_expr, na.rm = TRUE)
      
      data.frame(
        feature = unique(sub$feature),
        T_obs = T_obs,
        p_value = p_value,
        cct_value = cct_value,
        prob_sig = prob_sig,
        ave_expr = ave_expr,
        stringsAsFactors = FALSE
      )
    }
  )

  # Adjust p-values (BH)
  agg_df <- do.call(rbind, agg_list)
  agg_df$p_adj <- p.adjust(agg_df$p_value, method = "BH")
  agg_df$cct_adj <- p.adjust(agg_df$cct_value, method = "BH")
  
  # Final column order
  agg_df <- agg_df[
    c(
      "feature", 
      "T_obs", 
      "p_value",
      "p_adj",
      "cct_value",
      "cct_adj",
      "prob_sig", 
      "ave_expr"
    )
  ]
  rownames(agg_df) <- NULL
  
  return(agg_df)
}

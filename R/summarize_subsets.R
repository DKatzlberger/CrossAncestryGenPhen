#' Summarize subsets results by feature, with Cauchy-combined p-value
#'
#' @param data A data frame with columns: `feature`, `T_obs`, `p_value`, `p_adj`, `ave_expr`, `iteration`.
#'
#' @export
summarize_subsets <- function(
  data
) { 

  ## --- CTT helper function ---
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
  
  ## --- Aggregation ---
  feat_levels <- unique(data$feature)
  f <- factor(data$feature, levels = feat_levels, ordered = TRUE)

  parts <- split(data, f, drop = TRUE)

  agg_list <- lapply(parts, function(sub) {
    T_obs     <- mean(sub$T_obs,   na.rm = TRUE)
    p_value   <- mean(sub$p_value, na.rm = TRUE)
    cct_value <- cct_pval(sub$p_value, na.rm = TRUE)
    prob_sig  <- mean(sub$p_adj < 0.05, na.rm = TRUE)
    ave_expr  <- mean(sub$ave_expr, na.rm = TRUE)

    data.frame(
      feature   = as.character(sub$feature[1]),
      T_obs     = T_obs,
      p_value   = p_value,
      cct_value = cct_value,
      prob_sig  = prob_sig,
      ave_expr  = ave_expr,
      stringsAsFactors = FALSE
    )
  })

  agg_df <- do.call(rbind, agg_list)


  ## --- Adjust p-values (BH) ---
  agg_df$p_adj  <- p.adjust(agg_df$p_value,  method = "BH")
  agg_df$cct_adj<- p.adjust(agg_df$cct_value, method = "BH")
  

  ## --- Final column order ---
  agg_df <- agg_df[, c("feature","T_obs","p_value","p_adj",
                       "cct_value","cct_adj","prob_sig","ave_expr")]
  rownames(agg_df) <- NULL
  
  return(agg_df)
}

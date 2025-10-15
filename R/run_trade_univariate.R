#' Run univariate TRADE analysis for one contrast
#'
#' @param stats Data frame of stats with required cols.
#' @param coef_id Contrast ID to analyze (matches coef_id col).
#'
#' @return TRADE result object (list).
#'
#'
#'
#' @importFrom TRADEtools TRADE
#' @export
run_trade_univariate <- function(
  stats,
  coef_id = NULL
){

  ## --- Rename coef_id ---
  coef_id_ <- coef_id

  ## --- Check input ---
  required <- c(
    "coef_id", "coef_type", "contrast", "g_1", "g_2", "a_1", "a_2", 
    "feature", "T_obs", "SE", "p_value", "p_adj", "ave_expr"
  )

  # Stop if required columns are missing
  missing  <- setdiff(required, colnames(stats))
  if (length(missing) > 0) {
    stop("[run_trade_univariate] Missing required column(s): ", paste(missing, collapse = " "))
  }


  ## --- Available contrasts ---
  available <- unique(stats[, c("coef_id", "contrast")])

  if (!coef_id_ %in% stats$coef_id) {
    message("Contrast '", coef_id_, "' not found.\n")
    message("Available contrasts are:")

    lhs <- paste0(available$coef_id, ":")      
    width <- max(nchar(lhs))              

    for (i in seq_len(nrow(available))) {
      message(sprintf("%-*s %s", width, lhs[i], available$contrast[i]))
    }

    stop("[trade_univariate] Wrong contrast specified.", call. = FALSE)
  }


  ## --- Subset data for this contrast ---
  df <- stats[stats$coef_id == coef_id_, , drop = FALSE]


  ## --- Rename columns for TRADE ---
  trade_input <- data.frame(
    log2FoldChange = df$T_obs,
    lfcSE          = df$SE,
    pvalue         = df$p_value
  )
  rownames(trade_input) <- df$feature


  ## --- Run TRADE ---
  fit <- TRADEtools::TRADE(
    mode = "univariate",
    results1 = trade_input,
    annot_table = NULL,
    genes_exclude = NULL,
    n_sample = NULL
  )


  ## --- Return ---
  return(fit)
}

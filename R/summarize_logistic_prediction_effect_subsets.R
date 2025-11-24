#' Summarize logistic prediction effect subsets
#'
#' Computes per-iteration prediction performance (AUC or log-loss) for two
#' coefficient-defined subsets (`X` and `Y`) and summarizes their differences.
#'
#' @param stats A data frame containing at least the columns: `coef_id`, `coef_type`, `contrast`, `g_1`, `g_2`, `a_1`, `a_2`, `true`, `prob`, `iteration`.
#' @param method Performance metric to compute: `"auc"` or `"logloss"`.
#' @param by Optional character vector of additional grouping variables.
#'
#' @return A list with per-iteration correlations and summarized results.
#'
#' @importFrom data.table as.data.table
#' @importFrom data.table .SD
#' @importFrom yardstick roc_auc_vec mn_log_loss_vec
#' @importFrom stats quantile
#'
#' @export
summarize_logistic_prediction_effect_subsets <- function(
  stats,
  method = c("auc", "logloss"),
  by = NULL
){
  ## --- Match method ---
  method <- match.arg(method)


  ## --- Required columns ---
  required <- c("coef_id", "coef_type", "contrast", "g_1", "g_2", "a_1", "a_2", "true", "prob", "iteration")
  missing  <- setdiff(required, colnames(stats))
  if (length(missing) > 0) {
    stop("[summarize_logistic_prediction_effect_subsets] Missing required column(s): ", paste(missing, collapse = " "))
  }


  ## --- Check grouping vars ---
  if (!is.null(by) && !all(by %in% colnames(stats))) {
    stop("[summarize_logistic_prediction_effect_subsets] Grouping column(s) not found: ", paste(setdiff(by, colnames(stats)), collapse = " "))
  }

  # Group by iteration + user-defined grouping
  per_group <- unique(c("iteration", by))


  ## --- Metric functions ---
  methods_fun <- list(
    auc     = function(y, p) yardstick::roc_auc_vec(y, p, event_level = "second"),
    logloss = function(y, p) yardstick::mn_log_loss_vec(y, p, event_level = "second")
  )


  ## --- Data.table prep ---
  dt <- as.data.table(stats)
  # Metadata columns to keep 
  keep_cols <- c("coef_type", "g_1", "g_2", "a_1", "a_2")


  ## --- Compute per-iteration performance ---
  results <- lapply(
    names(methods_fun),
    function(m){
      ## Match the method
      fun <- methods_fun[[m]]

      ## Aggregate using fun
      methods_dt <- dt[, {
        cur <- .SD

        ## Subsets
        X <- cur[grepl("X$", coef_id)]
        Y <- cur[grepl("Y$", coef_id)]

        ## Metric
        metric_safe <- function(d) {
          d$true <- droplevels(d$true)
          fun(d$true, d$prob)
        }

        ## Compute metric
        RX <- metric_safe(X)
        RY <- metric_safe(Y)
        # Delta
        delta <- RX - RY


        ## Build result
        res <- list(
          RX    = RX,
          RY    = RY,
          delta = delta
        )

        ## Attach metadata
        for (col in keep_cols) res[[col]] <- cur[[col]][1]

        ## Return
        res
      }, by = per_group]

      ## Force output order
      out_cols <- c("coef_type", "contrast", "g_1", "g_2", "a_1", "a_2", "RX", "RY", "delta", per_group)
      out_cols <- out_cols[out_cols %in% colnames(methods_dt)]
      as.data.frame(methods_dt[, ..out_cols])
    }
  )
  names(results) <- names(methods_fun)

  ## --- Quantile CI ---
  summary_stats <- lapply(
    names(results), 
    function(m){
      tmp <- as.data.table(results[[m]])

      ## Compute summary stats per grouping combination
      summ <- tmp[, {
        d <- delta[!is.na(delta)]

        list(
          contrast   = unique(paste0(a_1, " - ", a_2)),
          delta_mean = mean(d),
          delta_q025 = quantile(d, 0.025),
          delta_q975 = quantile(d, 0.975)
        )
      }]

      ## Attach metadata columns (fixed across groups)
      for (col in keep_cols[keep_cols %in% colnames(dt)]) summ[[col]] <- dt[[col]][1]

      ## Force consistent output structure
      out_cols <- c("coef_type", "g_1", "g_2", "a_1", "a_2", "contrast", "delta_mean", "delta_q025", "delta_q975")
      out_cols <- out_cols[out_cols %in% colnames(summ)]
      as.data.frame(summ[, ..out_cols])
    }
  )
  names(summary_stats) <- names(methods_fun)

  ## --- Return ---
  list(
    all_res       = results,               
    all_delta_res = summary_stats,           
    sel_res       = results[[method]],    
    sel_delta_res = summary_stats[[method]]
  )
}
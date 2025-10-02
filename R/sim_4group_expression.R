#' Simulate RNA-seq Expression Data for a 2×2 Design
#'
#' Generate simulated RNA-seq count data for two ancestries, each with two
#' experimental groups. The function calls \code{\link{sim_2group_expression}}
#' for each ancestry and derives true baseline, relationship, and interaction
#' contrasts at the gene level.
#'
#' @param estimates_X Parameter estimates (means, dispersions, library sizes) for ancestry 1.
#' @param estimates_Y Parameter estimates for ancestry 2.
#' @param g_col Character scalar; name of the group column in metadata.
#' @param g_levels Character vector of length 2; group labels.
#' @param a_col Character scalar; name of the ancestry column in metadata.
#' @param a_levels Character vector of length 2; ancestry labels.
#' @param n_samples Integer; number of samples to simulate per condition.
#' @param n_degs Integer; number of DEGs per ancestry simulation.
#' @param log2fc Numeric; log2 fold-change magnitude for DEGs.
#' @param mean_method Method for extracting mean parameters.
#' @param disp_method Method for extracting dispersions.
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Logical; print a simulation summary.
#'
#'
#' @export
sim_4group_expression <- function(
  estimates_X,
  estimates_Y,
  g_col,
  g_levels,
  a_col,
  a_levels,
  n_samples,
  n_degs,
  log2fc,
  mean_method = c("mle", "map", "libnorm_mle", "libnorm_map"),
  disp_method = c("mle", "map"),
  seed = NULL,
  verbose = TRUE
){

  ## --- Match methods ---
  mean_method <- match.arg(mean_method)
  disp_method <- match.arg(disp_method)


  ## --- Factor levels ----
  if (length(g_levels) != 2 || length(a_levels) != 2) {
    stop("[sim_4group_expression] Function currently supports only 2x2 designs (two levels in g_level × two levels in a_level).")
  }

  a_1 <- a_levels[1]; a_2 <- a_levels[2]
  g_1 <- g_levels[1]; g_2 <- g_levels[2]


  ## --- Simulate two ancestries ----
  sim_X <- sim_2group_expression(
    estimates   = estimates_X,
    g_col       = g_col,
    g_levels    = g_levels,
    a_col       = a_col,
    a_level     = a_1,
    n_samples   = n_samples,
    n_degs      = n_degs,
    log2fc      = log2fc,
    mean_method = mean_method,
    disp_method = disp_method,
    seed        = seed
  )

  sim_Y <- sim_2group_expression(
    estimates   = estimates_Y,
    g_col       = g_col,
    g_levels    = g_levels,
    a_col       = a_col,
    a_level     = a_2,
    n_samples   = n_samples,
    n_degs      = n_degs,
    log2fc      = log2fc,
    mean_method = mean_method,
    disp_method = disp_method,
    seed        = if (is.null(seed)) NULL else seed + 1000
  )


  ## --- Extract truth ---
  # Per-condition means
  mean_X <- sim_X$feat
  mean_Y <- sim_Y$feat
  # Relationship (DEGs)
  rel_X <- sim_X$degs
  rel_Y <- sim_Y$degs
  # Features
  features <- mean_X$feature


  # Baseline_1: g1.a2 - g1.a1
  base_1 <- data.frame(
    coef_id   = "baseline_1",
    coef_type = "baseline",
    contrast  = paste0(g_1, ".", a_2, " - ", g_1, ".", a_1),
    g_1       = g_1, 
    g_2       = g_2,
    a_1       = a_1, 
    a_2       = a_2,
    feature   = features,
    T_obs     = log2(mean_Y$mean_g_1) - log2(mean_X$mean_g_1),
    ave_expr  = rowMeans(log2(cbind(mean_X$mean_g_1, mean_Y$mean_g_1) + 1)),
    is_DE     = as.numeric((log2(mean_Y$mean_g_1) - log2(mean_X$mean_g_1)) != 0),
    row.names = NULL
  )


  # Baseline_2: g2.a2 - g2.a1
  base_2 <- data.frame(
    coef_id   = "baseline_2",
    coef_type = "baseline",
    contrast  = paste0(g_2, ".", a_2, " - ", g_2, ".", a_1),
    g_1       = g_1, 
    g_2       = g_2,
    a_1       = a_1, 
    a_2       = a_2,
    feature   = features,
    T_obs     = log2(mean_Y$mean_g_2) - log2(mean_X$mean_g_2),
    ave_expr  = rowMeans(log2(cbind(mean_X$mean_g_2, mean_Y$mean_g_2) + 1)),
    is_DE     = as.numeric((log2(mean_Y$mean_g_2) - log2(mean_X$mean_g_2)) != 0),
    row.names = NULL
  )


  # Relationship_1: g2.a1 - g1.a1
  rel_1 <- rel_X
  rel_1$coef_id  <- "relationship_1"
  rel_1$coef_type <- "relationship"
  rel_1$a_1 <- a_1
  rel_1$a_2 <- a_1
  # Column order
  rel_1 <- rel_1[, c("coef_id","coef_type","contrast","g_1","g_2","a_1","a_2", "feature","T_obs","ave_expr","is_DE")]

  # Relationship_2: g2.a2 - g1.a2
  rel_2 <- rel_Y
  rel_2$coef_id  <- "relationship_2"
  rel_2$coef_type <- "relationship"
  rel_2$a_1 <- a_2
  rel_2$a_2 <- a_2
  # Column order
  rel_2 <- rel_2[, c("coef_id","coef_type","contrast","g_1","g_2","a_1","a_2", "feature","T_obs","ave_expr","is_DE")]


  # Interaction: (rel2 - rel1)
  int <- data.frame(
    coef_id   = "interaction",
    coef_type = "interaction",
    contrast  = paste0("(", g_2, ".", a_2, " - ", g_1, ".", a_2, ") - (", g_2, ".", a_1, " - ", g_1, ".", a_1, ")"),
    g_1       = g_1, 
    g_2       = g_2,
    a_1       = a_1, 
    a_2       = a_2,
    feature   = features,
    T_obs     = rel_2$T_obs - rel_1$T_obs,
    ave_expr  = rowMeans(log2(cbind(mean_X$mean_g_1, mean_X$mean_g_2, mean_Y$mean_g_1, mean_Y$mean_g_2) + 1)),
    is_DE     = as.numeric((rel_2$T_obs - rel_1$T_obs) != 0),
    row.names = NULL
  )


  ## --- Combine DEGs truths ---
  sim_degs <- rbind(base_1, base_2, rel_1, rel_2, int)


  ## --- Verbose message ---
  if (verbose) {
    fmt_counts <- function(meta, label) {
      tab <- table(meta[[g_col]], dnn = NULL)
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = " ")
    }

    message("\n4-group simulation summary:")
    message(sprintf("%s (X):    N: %-4d  n_DEGs: %-4d  log2FC: %-4.1f %s features: %-4d", a_1, nrow(sim_X$matr), n_degs, log2fc, fmt_counts(sim_X$meta, a_1), ncol(sim_X$matr)))
    message(sprintf("%s (Y):    N: %-4d  n_DEGs: %-4d  log2FC: %-4.1f %s features: %-4d", a_2, nrow(sim_Y$matr), n_degs, log2fc, fmt_counts(sim_Y$meta, a_2), ncol(sim_X$matr)))
  }


  ## --- Return ---
  return(
    list(
      X = list(
        matr = sim_X$matr, 
        meta = sim_X$meta
        # feat = sim_X$feat
        # plot = sim_X$in_out_plots
      ),
      Y = list(
        matr = sim_Y$matr, 
        meta = sim_Y$meta
        # feat = sim_Y$feat
        # plot  = sim_Y$in_out_plots
      ),
      degs = sim_degs
    )
  )
}
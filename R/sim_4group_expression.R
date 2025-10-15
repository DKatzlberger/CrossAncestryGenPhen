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
  mean_method = c("mle", "libnorm_mle"),
  disp_method = c("mle"),
  drop_zeros = FALSE,
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


  ## --- Drop problematic features ---
  if (drop_zeros){
    # Pick means & dispersions consistently
    means_X <- estimates_X$means[[mean_method]]
    means_Y <- estimates_Y$means[[mean_method]]
    disps_X <- estimates_X$disps[[disp_method]]
    disps_Y <- estimates_Y$disps[[disp_method]]

    bad_idx <- which(means_X == 0 | means_Y == 0 | disps_X == 0 | disps_Y == 0)

    if (length(bad_idx) > 0) {
      message(sprintf("\n[sim_4group_expression] Dropping %d feature(s) with zero mean/disp before simulation.", length(bad_idx)))

      keep <- setdiff(seq_along(means_X), bad_idx)

      # Subset both estimate objects
      estimates_X$means <- lapply(estimates_X$means, `[`, keep)
      estimates_X$disps <- lapply(estimates_X$disps, `[`, keep)

      estimates_Y$means <- lapply(estimates_Y$means, `[`, keep)
      estimates_Y$disps <- lapply(estimates_Y$disps, `[`, keep)
    }
  }


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

  # Feature alignment
  if (!identical(mean_X$feature, mean_Y$feature)) {
    stop("[sim_4group_expression] Feature means are not identical for X and Y.")
  }

  # Relationship (DEGs)
  rel_X <- sim_X$degs
  rel_Y <- sim_Y$degs

  # Feature alignment
  if (!identical(rel_X$feature, rel_Y$feature)) {
    stop("[sim_4group_expression] DEGs are not identical for X and Y.")
  }

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
    row.names = NULL
  )


  # Relationship_1: g2.a1 - g1.a1
  rel_X$coef_id   <- "relationship_1"
  rel_X$coef_type <- "relationship"
  rel_X$a_1       <- a_1
  rel_X$a_2       <- a_1
  # Column order
  rel_X <- rel_X[, c("coef_id", "coef_type", "contrast", "g_1", "g_2", "a_1", "a_2", "feature", "T_obs")]

  # Relationship_2: g2.a2 - g1.a2
  rel_Y$coef_id   <- "relationship_2"
  rel_Y$coef_type <- "relationship"
  rel_Y$a_1       <- a_2
  rel_Y$a_2       <- a_2
  # Column order
  rel_Y <- rel_Y[, c("coef_id", "coef_type", "contrast", "g_1", "g_2", "a_1", "a_2", "feature", "T_obs")]


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
    T_obs     = rel_Y$T_obs - rel_X$T_obs,
    row.names = NULL
  )


  ## --- Combine DEGs truths ---
  sim_degs <- rbind(base_1, base_2, rel_X, rel_Y, int)


  ## --- Sanity check: look for NA in truth table ---
  if (anyNA(sim_degs)) {
    bad_rows <- which(!complete.cases(sim_degs))
    warning(sprintf("\n[sim_4group_expression] %d row(s) in truth table contain NA values.", length(bad_rows)))
  }


  ## --- Align to DEG contrasts ---
  {
    # Build "true" group means table
    M <- data.frame(
      g1.a1 = mean_X$mean_g_1,
      g2.a1 = mean_X$mean_g_2,
      g1.a2 = mean_Y$mean_g_1,
      g2.a2 = mean_Y$mean_g_2
    )

    # Compute DGE contrasts directly from means
    dge_truth <- data.frame(
      baseline_1     = log2(M$g1.a2) - log2(M$g1.a1),
      baseline_2     = log2(M$g2.a2) - log2(M$g2.a1),
      relationship_1 = log2(M$g2.a1) - log2(M$g1.a1),
      relationship_2 = log2(M$g2.a2) - log2(M$g1.a2),
      interaction    = (log2(M$g2.a2) - log2(M$g1.a2)) - (log2(M$g2.a1) - log2(M$g1.a1))
    )

    # Corresponding simulation truth
    sim_truth <- data.frame(
      baseline_1     = log2(mean_Y$mean_g_1) - log2(mean_X$mean_g_1),
      baseline_2     = log2(mean_Y$mean_g_2) - log2(mean_X$mean_g_2),
      relationship_1 = rel_X$T_obs,
      relationship_2 = rel_Y$T_obs,
      interaction    = rel_Y$T_obs - rel_X$T_obs
    )

    # Check numeric alignment
    diffs <- dge_truth - sim_truth
    max_diff <- sapply(diffs, function(x) max(abs(x), na.rm = TRUE))
    corrs <- mapply(function(x, y) suppressWarnings(cor(x, y, use = "pairwise.complete.obs")), dge_truth, sim_truth)

    bad <- which(max_diff > 1e-10 | corrs < 0.999)
    if (length(bad) > 0) {
      stop(sprintf("[sim_4group_expression] Contrast alignment failed for: %s", paste(names(bad), collapse = ", ")))
    }
  }


  ## --- Final consistency checks ----
  {
    # Check sample alignment
    if (!identical(rownames(sim_X$matr), rownames(sim_X$meta)) ||
        !identical(rownames(sim_Y$matr), rownames(sim_Y$meta))) {
      stop("[sim_4group_expression] Sample IDs in 'matrix' and 'meta' not aligned.")
    }

    # Check feature alignment
    feat_matr_X <- colnames(sim_X$matr)
    feat_matr_Y <- colnames(sim_Y$matr)
    feat_truth  <- unique(sim_degs$feature)

    if (!identical(feat_matr_X, feat_matr_Y) || !identical(feat_matr_X, feat_truth)) {
      stop("[sim_4group_expression] Feature names misaligned between 'matrix' and truth tables.")
    }
  }


  ## --- Plot ---



  ## --- Verbose message ---
  if (verbose) {

    fmt_counts <- function(meta, g_col) {
      tab <- table(meta[[g_col]], dnn = NULL)
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = "  ")
    }

    # Define combined groups
    grp_X <- interaction(sim_X$meta[[g_col]], sim_X$meta[[a_col]], drop = TRUE)
    grp_Y <- interaction(sim_Y$meta[[g_col]], sim_Y$meta[[a_col]], drop = TRUE)
    grp   <- c(grp_X, grp_Y)

    message("\n4-group simulation summary:")
    message(sprintf("%-18s %s", "Groups:", paste(unique(grp), collapse = "  ")))

    message(sprintf(
      "%-18s N: %-5d  n_DEGs: %-18d  log2FC: %-18.1f  %-18s  features: %-18d",
      paste0(a_1, " (X):"),
      nrow(sim_X$matr),
      n_degs,
      log2fc,
      fmt_counts(sim_X$meta, g_col),
      ncol(sim_X$matr)
    ))

    message(sprintf(
      "%-18s N: %-18d  n_DEGs: %-18d  log2FC: %-18.1f  %-18s  features: %-18d",
      paste0(a_2, " (Y):"),
      nrow(sim_Y$matr),
      n_degs,
      log2fc,
      fmt_counts(sim_Y$meta, g_col),
      ncol(sim_Y$matr)
    ))
  }


  # if (verbose) {
  #   fmt_counts <- function(meta, label) {
  #     tab <- table(meta[[g_col]], dnn = NULL)
  #     paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = " ")
  #   }

  #   # Define group labels for clarity
  #   grp_X <- interaction(sim_X$meta[[g_col]], sim_X$meta[[a_col]], drop = TRUE)
  #   grp_Y <- interaction(sim_Y$meta[[g_col]], sim_Y$meta[[a_col]], drop = TRUE)
  #   grp   <- c(grp_X, grp_Y)

  #   message("\n4-group simulation summary:")
  #   message(sprintf("Groups:    %s", paste(unique(grp), collapse = "  ")))
  #   message(sprintf("%s (X):    N: %-4d  n_DEGs: %-4d  log2FC: %-4.1f %s features: %-4d", a_1, nrow(sim_X$matr), n_degs, log2fc, fmt_counts(sim_X$meta, a_1), ncol(sim_X$matr)))
  #   message(sprintf("%s (Y):    N: %-4d  n_DEGs: %-4d  log2FC: %-4.1f %s features: %-4d", a_2, nrow(sim_Y$matr), n_degs, log2fc, fmt_counts(sim_Y$meta, a_2), ncol(sim_X$matr)))
  # }


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
        # feat  = sim_Y$feat
        # plot  = sim_Y$in_out_plots
      ),
      degs = sim_degs
    )
  )
}
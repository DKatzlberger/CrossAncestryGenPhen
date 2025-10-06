#' Simulate RNA-seq Expression Data for Two Groups
#'
#' Generate a simulated RNA-seq count matrix for two groups using parameters
#' estimated from real data and \code{compcodeR::generateSyntheticData}.
#' Simulated data retain similar mean–variance characteristics as the input,
#' with a specified number of differentially expressed genes (DEGs).
#'
#' 
#' @param estimates Optional list of pre-computed NB parameter estimates.
#' @param g_col Character scalar giving the name of the group column in the sample metadata.
#' @param g_levels Character vector of length 2 giving the group labels.
#' @param a_col Character scalar giving the name of the ancestry column in the sample metadata.
#' @param a_level Character scalar giving the ancestry label for this simulation.
#' @param n_samples Integer, number of samples to simulate per condition.
#' @param n_degs Integer, number of differentially expressed genes to simulate.
#' @param log2fc Numeric, log2 fold-change magnitude for DEGs.
#' @param mean_method Character string, method to use for mean estimates.
#' @param disp_method Character string, method to use for dispersion estimates.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @export
#'
#' @importFrom compcodeR generateSyntheticData
#' @importFrom stats setNames
sim_2group_expression <- function(
  estimates,
  g_col,
  g_levels,
  a_col,
  a_level,
  n_samples,
  n_degs,
  log2fc,
  mean_method = c("mle", "map", "libnorm_mle", "libnorm_map"),
  disp_method = c("mle", "map"),
  seed = NULL 
){

  ## --- Set seed ---
  if (!is.null(seed)) set.seed(seed)


  ## --- Match method arguments safely ---
  mean_method <- match.arg(mean_method)
  disp_method <- match.arg(disp_method)

  ## --- Phenotype factors ----
  if (length(g_levels) != 2 || length(a_level) != 1) {
    stop("[sim_2group_expression] Function currently supports only 2x1 designs (two levels in g_level × one level in a_level).")
  }
  g_1 <- g_levels[1]; g_2 <- g_levels[2]


  ## ---  Estimate parameters ---
  real_estimates <- estimates


  ## --- Extract means ---
  real_means <- switch(
    mean_method,
    mle = real_estimates$means$mle,
    map = real_estimates$means$map,
    libnorm_mle = real_estimates$means$libnorm_mle,
    libnorm_map = real_estimates$means$libnorm_map
  )


  ## --- Extract dispersions ---
  real_disps <- switch(
    disp_method,
    mle = real_estimates$disps$mle,
    map = real_estimates$disps$map
  )


  ## --- Check for zero means ---
  n_zero_means <- sum(real_means == 0 | is.na(real_means))
  if (n_zero_means > 0) {
    warning(sprintf("\n[sim_2group_expression] %d feature(s) with zero mean detected; these may cause NaN/NA logFC or is_DE values in truth table.", n_zero_means))
  }

  ## --- Check for zero dispersions ---
  n_zero_disps <- sum(real_disps == 0 | is.na(real_disps))
  if (n_zero_disps > 0) {
    warning(sprintf("\n[sim_2group_expression] %d feature(s) with zero dispersion detected; simulation may generate degenerate counts.", n_zero_disps))
  }


  ## --- Use generateSyntheticData from compcodeR ---
  sim <- compcodeR::generateSyntheticData(
    dataset = a_level,
    samples.per.cond = n_samples,
    relmeans = real_means,
    dispersions = real_disps,
    n.vars = length(real_means),
    n.diffexp = n_degs,
    effect.size = if (log2fc == 0) 0 else 2^log2fc,
    fraction.upregulated = 0.5,
    seqdepth = real_estimates$libsize$mean,
    minfact  = real_estimates$libsize$min / real_estimates$libsize$mean,
    maxfact  = real_estimates$libsize$max / real_estimates$libsize$mean,
    fraction.non.overdispersed = 0,
    random.outlier.high.prob = 0,
    random.outlier.low.prob = 0,
    single.outlier.high.prob = 0,
    single.outlier.low.prob = 0,
    filter.threshold.total = 0,
    filter.threshold.mediancpm = 0,
    repl.id = seed
  )


  ## --- Extract counts (sample x features) ---
  sim_counts <- sim@count.matrix
  colnames(sim_counts) <- paste0(a_level, "_", seq_len(ncol(sim_counts)))
  sim_counts <- t(sim_counts) # rownames are samples (matrix)


  ## --- Meta ---
  sim_meta <- sim@sample.annotations
  sim_meta[[g_col]]  <- factor(sim_meta$condition, levels = c(1, 2), labels = c(g_1, g_2))
  sim_meta[[a_col]]  <- a_level
  rownames(sim_meta) <- rownames(sim_counts) # rownames are samples (matrix)

  
  ## --- Group means ---
  sim_mean <- data.frame(
    g_1       = g_1,
    g_2       = g_2,
    mean_g_1  = sim@variable.annotations$truemeans.S1,
    mean_g_2  = sim@variable.annotations$truemeans.S2,
    feature   = colnames(sim_counts),
    row.names = NULL
  )


  ## --- Relationship ---
  sim_degs <- data.frame(
    contrast  = paste0(g_2, ".", a_level, " - ", g_1, ".", a_level),
    g_1       = g_1,
    g_2       = g_2,
    feature   = colnames(sim_counts),
    T_obs     = sim@variable.annotations$truelog2foldchanges,
    ave_expr  = rowMeans(log2(cbind(sim_mean$mean_g_1, sim_mean$mean_g_2) + 1)),
    is_DE     = sim@variable.annotations$differential.expression, 
    row.names = NULL
  )
  # rownames(sim_features) <- colnames(sim_counts) # rownames are feature names (matrix)


  # ## --- Comparison real vs. sim ---
  # sim_estimates <- estimate_params(sim_counts, verbose = FALSE)


  # ## --- Gene means ---
  # real_sim_means <- plot_estimated_means(
  #   estimates_X = real_estimates,
  #   estimates_Y = sim_estimates,
  #   method = "mle",
  #   ancestry_X = a_level, 
  #   ancestry_Y = sim_ancestry,
  #   title = "Simulated vs real gene-wise means"
  # )

  # ## --- Gene dispersions ---
  # real_sim_disps <- plot_estimated_dispersions(
  #   estimates_X = real_estimates,
  #   estimates_Y = sim_estimates,
  #   method = "mle",
  #   ancestry_X = a_level, 
  #   ancestry_Y = sim_ancestry,
  #    title = "Simulated vs real gene-wise dispersions"
  # )


  # ## --- Plots ---
  # plots <- list(
  #   means = real_sim_means,
  #   disps = real_sim_disps
  # )


  ## --- Return ---
  return(
    list(
      matr = sim_counts,
      meta = sim_meta,
      feat = sim_mean,
      degs = sim_degs
      # input_params  = real_estimates,
      # output_params = sim_estimates,
      # in_out_plots  = plots
    )
  )
}
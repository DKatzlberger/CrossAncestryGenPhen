#' Simulate RNA-seq Expression Data for Two Groups
#'
#' Generate a simulated RNA-seq count matrix for two groups using parameters
#' estimated from real data and \code{compcodeR::generateSyntheticData}.
#' Simulated data retain similar mean–variance characteristics as the input,
#' with a specified number of differentially expressed genes (DEGs).
#'
#' @param X Numeric matrix or data frame of counts from the real data
#'        (samples in rows, genes in columns).
#' @param estimates Optional list of pre-computed parameter estimates
#'        from \code{\link{estimate_params}}. If \code{NULL}, parameters
#'        are estimated from \code{X}.
#' @param ancestry Character scalar giving the ancestry label for this simulation.
#' @param n_samples Integer, number of samples to simulate per condition.
#' @param n_degs Integer, number of differentially expressed genes to simulate.
#' @param log2FC Numeric, log2 fold-change magnitude for DEGs.
#' @param mean_method Character string, method to use for mean estimates.
#'        One of \code{"mle"}, \code{"map"}, \code{"libnorm_mle"}, \code{"libnorm_map"}.
#' @param disp_method Character string, method to use for dispersion estimates.
#'        One of \code{"mle"}, \code{"map"}.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{X}}{Simulated count matrix (samples x genes).}
#'   \item{\code{M}}{Data frame of sample metadata.}
#'   \item{\code{f}}{Data frame of gene-level features:
#'        DE status (\code{is_DE}) and true log2FC (\code{true_log2FC}).}
#'   \item{\code{input_params}}{Parameter estimates from the real data.}
#'   \item{\code{output_params}}{Parameter estimates from the simulated data.}
#'   \item{\code{in_out_plots}}{List of ggplot objects comparing means and dispersions
#'        between real and simulated data.}
#' }
#' @export
#'
#' @importFrom compcodeR generateSyntheticData
#' @importFrom stats setNames
sim_2group_expression <- function(
  estimates,
  ancestry,
  n_samples,
  n_degs,
  log2FC,
  mean_method = c("mle", "map", "libnorm_mle", "libnorm_map"),
  disp_method = c("mle", "map"),
  seed = NULL 
){

  ## --- Set seed ---
  if (!is.null(seed)) set.seed(seed)


  ## --- Match method arguments safely ---
  mean_method <- match.arg(mean_method)
  disp_method <- match.arg(disp_method)


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


  ## --- Use generateSyntheticData from compcodeR ---
  sim_ancestry <- paste0(ancestry, "_sim")
  sim <- compcodeR::generateSyntheticData(
    dataset = sim_ancestry,
    samples.per.cond = n_samples,
    relmeans = real_means,
    dispersions = real_disps,
    n.vars = length(real_means),
    n.diffexp = n_degs,
    effect.size = if (log2FC == 0) 0 else 2^log2FC,
    fraction.upregulated = 0.5,
    seqdepth = real_estimates$libsize$mean,
    minfact = real_estimates$libsize$min,
    maxfact = real_estimates$libsize$max,
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
  colnames(sim_counts) <- paste0(sim_ancestry, "_", seq_len(ncol(sim_counts)))
  sim_counts <- t(sim_counts) # rownames are samples (matrix)


  ## --- Meta ---
  sim_meta <- sim@sample.annotations
  sim_meta$condition <- factor(sim_meta$condition, levels = c(1, 2), labels = c("Control", "Case"))
  sim_meta$ancestry <- sim_ancestry
  rownames(sim_meta) <- rownames(sim_counts) # rownames are samples (matrix)


  ## --- Features ---
  is_DE <- sim@variable.annotations$differential.expression
  true_log2FC <- sim@variable.annotations$truelog2foldchanges
  sim_features <- data.frame(is_DE = is_DE, true_log2FC = true_log2FC)
  rownames(sim_features) <- colnames(sim_counts) # rownames are feature names (matrix)


  ## --- Comparison real vs. sim ---
  sim_estimates <- estimate_params(sim_counts)


  ## --- Gene means ---
  real_sim_means <- plot_estimated_means(
    estimates_X = real_estimates,
    estimates_Y = sim_estimates,
    method = "mle",
    ancestry_X = ancestry, 
    ancestry_Y = sim_ancestry,
    title = "Simulated vs real gene-wise means"
  )

  ## --- Gene dispersions ---
  real_sim_disps <- plot_estimated_dispersions(
    estimates_X = real_estimates,
    estimates_Y = sim_estimates,
    method = "mle",
    ancestry_X = ancestry, 
    ancestry_Y = sim_ancestry,
     title = "Simulated vs real gene-wise dispersions"
  )


  ## --- Plots ---
  plots <- list(
    means = real_sim_means,
    disps = real_sim_disps
  )


  ## --- Return ---
  return(
    list(
      counts        = sim_counts,
      meta          = sim_meta,
      feats         = sim_features,
      input_params  = real_estimates,
      output_params = sim_estimates,
      in_out_plots  = plots
    )
  )
}
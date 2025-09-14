#' Estimate RNA-seq parameters with edgeR
#'
#' Fits an intercept-only negative binomial GLM to RNA-seq count data using
#' \pkg{edgeR}. Returns mean expression, dispersion estimates, library sizes,
#' and optional QC plots.
#'
#' @param X Numeric matrix of raw counts (samples x genes).
#' @param plot Logical, default `TRUE`. If `TRUE`, print QC plots of mean
#'   expression, dispersion distribution, and the mean–dispersion trend.
#' @param seed Optional integer for reproducibility.
#' @param verbose Logical, default `TRUE`. Print a brief summary.
#'
#' @return A list with:
#' \describe{
#'   \item{mains}{Dataset info (samples, features, gene IDs).}
#'   \item{means}{Mean expression estimates (raw, logCPM, MLE, MAP, lib-normalized).}
#'   \item{disps}{Dispersion estimates (common, trended, MLE, MAP).}
#'   \item{libsize}{Effective library sizes (mean, min, max).}
#'   \item{plot}{Patchwork object with QC plots.}
#' }
#'
#' @details Library sizes are normalized with
#' \code{\link[edgeR]{calcNormFactors}} before dispersion estimation.
#'
#' @importFrom edgeR DGEList calcNormFactors estimateGLMCommonDisp
#'   estimateGLMTrendedDisp estimateGLMTagwiseDisp glmFit cpm
#' @import ggplot2
#' @import patchwork
#' @export
estimate_params <- function(
  X,
  plot = TRUE,
  seed = NULL,
  verbose = TRUE
) {

  ## --- Input data structure check ---
  assert_input(X = X)


  ## --- Seed ---
  if (!is.null(seed)) set.seed(seed)
  

  ## --- Transpose to genes x samples ---
  X <- t(X)


  ## --- Summary info ---
  mains <- list(
    n_samples = ncol(X),
    n_features = nrow(X),
    features = rownames(X)
  )


  ## --- Intercept-only model ---
  design <- matrix(1, ncol = 1, nrow = ncol(X))


  ## --- Create DGEList ---
  dge <- edgeR::DGEList(counts = X)


  ## --- Estimate library sizes ---
  dge <- edgeR::calcNormFactors(dge)
  eff_libsizes <- dge$samples$lib.size * dge$samples$norm.factors
  mean_eff_libsize <- mean(eff_libsizes)
  min_eff_libsize  <- min(eff_libsizes)
  max_eff_libsize  <- max(eff_libsizes)

  estimated_libsize <- list(
    mean = mean_eff_libsize,
    min = min_eff_libsize,
    max = max_eff_libsize
  )


  ## --- Estimate dispersions ---
  dge_disp <- edgeR::estimateGLMCommonDisp(dge, design = design)
  dge_disp <- edgeR::estimateGLMTrendedDisp(dge_disp, design = design)

  # MLE dispersion (no prior)
  dge_mle <- edgeR::estimateGLMTagwiseDisp(dge_disp, design = design, prior.df = 0)
  fit_mle <- edgeR::glmFit(dge_mle, design, dispersion = dge_mle$tagwise.dispersion)
  means_mle <- rowMeans(fit_mle$fitted.values)

  # MAP dispersion (with prior)
  dge_map <- edgeR::estimateGLMTagwiseDisp(dge_disp, design = design)
  fit_map <- edgeR::glmFit(dge_map, design, dispersion = dge_map$tagwise.dispersion)
  means_map <- rowMeans(fit_map$fitted.values)


  ## --- Estimate means ---
  # Raw and logCPM means
  means_raw <- rowMeans(X)
  means_logcpm <- rowMeans(edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE))

  # Lambda estimates (normalized by effective lib size)
  lambda_mle <- rowMeans(sweep(fit_mle$fitted.values, 2, eff_libsizes, FUN = "/"))
  lambda_map <- rowMeans(sweep(fit_map$fitted.values, 2, eff_libsizes, FUN = "/"))

  estimated_means <- list(
    raw = as.numeric(means_raw),
    logcpm = as.numeric(means_logcpm),
    mle = as.numeric(means_mle),
    map = as.numeric(means_map),
    libnorm_mle = as.numeric(lambda_mle),
    libnorm_map = as.numeric(lambda_map)
  )

  estimated_disps <- list(
    common = dge_disp$common.dispersion,
    trend = dge_disp$trended.dispersion,
    mle = dge_mle$tagwise.dispersion,
    map = dge_map$tagwise.dispersion
  )


  ## --- Verbose massage ---
  if (verbose){
    message("\nEstimate NB params summary:")
    message(sprintf("Dataset:        groups:  %-4d  N: %-4d  features: %-4d", ncol(design), ncol(X), nrow(X)))
    message(sprintf("Means (logCPM): mean:  %-4.1f  median: %-4.1f sd: %-4.1f", mean(means_logcpm), median(means_logcpm), sd(means_logcpm)))
    message(sprintf("Disps (MLE):    mean:  %-4.1f  median: %-4.1f sd: %-4.1f", mean(dge_mle$tagwise.dispersion), median(dge_mle$tagwise.dispersion), sd(dge_mle$tagwise.dispersion)))
    message(sprintf("Lib. size:      mean:  %-4.1e  min: %-4.1e  max: %-4.1e", mean_eff_libsize, min_eff_libsize, max_eff_libsize))
  }


  ## --- Plot ---
  data <- data.frame(
    means_logcpm = means_logcpm,
    disp_mle = dge_mle$tagwise.dispersion,
    bcv = sqrt(dge_mle$tagwise.dispersion)
  )

  # Mean distribution
  p_means <- ggplot(data, aes(x = means_logcpm)) +
    geom_histogram(bins = 30, fill = "grey80", color = "black", linewidth = 0.1) +
    geom_vline(aes(xintercept = mean(means_logcpm)), color = "red", linewidth = 0.3) +
    labs(x = "Log2 CPM", y = "Count")+ theme_nature_fonts() +
    theme_small_legend() + theme_white_background()

  # Dispersion distribution (MLE)
  p_disps <- ggplot(data, aes(x = log2(disp_mle))) +
    geom_histogram(bins = 30, fill = "grey80", color = "black", linewidth = 0.1) +
    geom_vline(aes(xintercept = mean(log2(disp_mle))), color = "red", linewidth = 0.3) +
    labs(x = "Log2 Dispersion", y = "Count") + theme_nature_fonts() +
    theme_small_legend() + theme_white_background()


  # Mean–dispersion relationship
  p_trend <- ggplot(data, aes(x = means_logcpm, y = log2(disp_mle))) +
    geom_point(size = 0.5) + geom_smooth(se = FALSE, color = "red", linewidth = 0.3) +
    labs(x = "Log2 CPM", y = "Log2 Dispersion") + theme_nature_fonts() +
    theme_small_legend() + theme_white_background()

  # Patchwork
  p <- patchwork::wrap_plots(p_means, p_disps, p_trend, ncol = 3, nrow = 1)
  if (plot){
    print(p)
  }



  ## --- Return ---
  return(
      list(
      mains = mains,
      means = estimated_means,
      disps = estimated_disps,
      libsize = estimated_libsize,
      plot = p
    )
  )
}
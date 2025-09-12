#' Estimate RNA-seq Model Parameters from Count Data
#'
#' Estimate gene-level mean expression and dispersion parameters from RNA-seq
#' count data using \pkg{edgeR}. The function fits an intercept-only GLM model
#' to compute maximum likelihood estimates (MLE) and maximum a posteriori
#' (MAP) estimates for means and dispersions, along with raw and logCPM means.
#'
#' @param X Numeric matrix of raw RNA-seq counts with \strong{samples in rows}
#'          and \strong{genes in columns}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{mains}}{List with basic dataset information:
#'     \code{n_samples} (integer), \code{n_features} (integer), and
#'     \code{features} (gene IDs).}
#'   \item{\code{means}}{List of mean expression estimates:
#'     \code{raw} (raw mean counts),
#'     \code{logcpm} (mean log2 CPM),
#'     \code{mle} (fitted values from MLE dispersion),
#'     \code{map} (fitted values from MAP dispersion),
#'     \code{libnorm_mle} (MLE fitted means normalized by effective library size),
#'     \code{libnorm_map} (MAP fitted means normalized by effective library size).}
#'   \item{\code{disps}}{List of dispersion estimates:
#'     \code{common} (common dispersion),
#'     \code{trend} (trended dispersion),
#'     \code{mle} (tagwise dispersion without prior),
#'     \code{map} (tagwise dispersion with prior).}
#'   \item{\code{libsize}}{Numeric scalar giving the mean effective library size.}
#' }
#'
#' @details
#' The function uses an intercept-only design matrix to estimate baseline
#' mean expression and dispersion parameters across all samples. The effective
#' library size is computed as the product of the raw library size and the
#' normalization factor estimated by \code{\link[edgeR]{calcNormFactors}}.
#'
#' \strong{Important:} The input \code{X} must be a matrix with samples in rows
#' and genes in columns. Internally, the function transposes \code{X} to match
#' the gene-by-sample format expected by \pkg{edgeR}.
#'
#' @seealso
#' \code{\link[edgeR]{DGEList}},
#' \code{\link[edgeR]{estimateGLMCommonDisp}},
#' \code{\link[edgeR]{estimateGLMTagwiseDisp}},
#' \code{\link[edgeR]{glmFit}}
#'
#'
#' @importFrom edgeR DGEList calcNormFactors estimateGLMCommonDisp
#'             estimateGLMTrendedDisp estimateGLMTagwiseDisp glmFit cpm
#'
#' @export
estimate_params <- function(
  X,
  seed = NULL
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

  return(list(
    mains = mains,
    means = estimated_means,
    disps = estimated_disps,
    libsize = estimated_libsize
  ))
}
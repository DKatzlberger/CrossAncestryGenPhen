#' Simulate RNA-seq Expression Data for Four Groups (Two Ancestries × Two Conditions)
#'
#' Generate simulated RNA-seq count matrices for two ancestries, each with two
#' conditions, based on parameters estimated from real data. This function calls
#' \code{\link{sim_2group_expression}} twice (once per ancestry) and then calculates
#' the true interaction effects between ancestries for each gene.
#'
#' @param X Numeric matrix or data frame of counts for the first ancestry
#'        (samples in rows, genes in columns).
#' @param Y Numeric matrix or data frame of counts for the second ancestry
#'        (samples in rows, genes in columns).
#' @param estimates_X Optional list of parameter estimates from \code{\link{estimate_params}}
#'        for ancestry X. If \code{NULL}, parameters are estimated from \code{X}.
#' @param estimates_Y Optional list of parameter estimates from \code{\link{estimate_params}}
#'        for ancestry Y. If \code{NULL}, parameters are estimated from \code{Y}.
#' @param ancestry_X Character scalar giving the ancestry label for \code{X}.
#' @param ancestry_Y Character scalar giving the ancestry label for \code{Y}.
#' @param n_samples_X Integer, number of samples to simulate per condition for ancestry X.
#' @param n_samples_Y Integer, number of samples to simulate per condition for ancestry Y.
#' @param n_degs_X Integer, number of differentially expressed genes to simulate in ancestry X.
#' @param n_degs_Y Integer, number of differentially expressed genes to simulate in ancestry Y.
#' @param log2FC_X Numeric, log2 fold-change magnitude for DEGs in ancestry X.
#' @param log2FC_Y Numeric, log2 fold-change magnitude for DEGs in ancestry Y.
#' @param mean_method Character string, method to use for mean estimates in both ancestries.
#'        One of \code{"mle"}, \code{"map"}, \code{"libnorm_mle"}, \code{"libnorm_map"}.
#' @param disp_method Character string, method to use for dispersion estimates in both ancestries.
#'        One of \code{"mle"}, \code{"map"}.
#' @param seed Optional integer random seed for reproducibility. The simulation
#'        for ancestry Y will use \code{seed + 1000} to ensure different DE gene sets.
#' @param verbose Logical; if `TRUE`, print a summary of simulated groups.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{X}}{Simulated count matrix (samples x genes) for ancestry X.}
#'   \item{\code{Y}}{Simulated count matrix (samples x genes) for ancestry Y.}
#'   \item{\code{MX}}{Sample metadata for ancestry X.}
#'   \item{\code{MY}}{Sample metadata for ancestry Y.}
#'   \item{\code{fX}}{Gene-level features (DE status, true log2FC) for ancestry X.}
#'   \item{\code{fY}}{Gene-level features (DE status, true log2FC) for ancestry Y.}
#'   \item{\code{pX}}{List of ggplot objects comparing means and dispersions for ancestry X.}
#'   \item{\code{pY}}{List of ggplot objects comparing means and dispersions for ancestry Y.}
#'   \item{\code{fI}}{Data frame of interaction effects, with DE status and true interaction log2FC.}
#' }
#'
#' @details
#' The interaction effect for each gene is defined as:
#' \deqn{\mathrm{Interaction\ log2FC} = \mathrm{log2FC}_Y - \mathrm{log2FC}_X}
#' where \eqn{\mathrm{log2FC}_X} and \eqn{\mathrm{log2FC}_Y} are the true log2
#' fold-changes from the simulated data for ancestry X and Y respectively.
#'
#' @seealso \code{\link{sim_2group_expression}}
#'
#' @export
sim_4group_expression <- function(
  estimates_X,
  estimates_Y,
  g_col,
  g_levels,
  a_col,
  a_levels,
  n_samples_X,
  n_samples_Y,
  n_degs_X,
  n_degs_Y,
  log2FC_X,
  log2FC_Y,
  mean_method = c("mle", "map", "libnorm_mle", "libnorm_map"),
  disp_method = c("mle", "map"),
  seed = NULL,
  verbose = TRUE
){

  ## --- Match methods ---
  mean_method <- match.arg(mean_method)
  disp_method <- match.arg(disp_method)


  ## --- Ancestry levels ----
  if (length(g_levels) != 2 || length(a_levels) != 2) {
    stop("[sim_4group_expression] Function currently supports only 2x2 designs (two levels in g_level × two levels in a_level).")
  }

  a_1 <- a_levels[1]; a_2 <- a_levels[2]
  a_1_sim <- paste0(a_1, "_sim"); a_2_sim <- paste0(a_2, "_sim")


  ## --- Simulate two ancestries ----
  sim_X <- sim_2group_expression(
    estimates   = estimates_X,
    g_col       = g_col,
    g_levels    = g_levels,
    a_col       = a_col,
    a_level     = a_1_sim,
    n_samples   = n_samples_X,
    n_degs      = n_degs_X,
    log2FC      = log2FC_X,
    mean_method = mean_method,
    disp_method = disp_method,
    seed        = seed
  )

  sim_Y <- sim_2group_expression(
    estimates   = estimates_Y,
    g_col       = g_col,
    g_levels    = g_levels,
    a_col       = a_col,
    a_level     = a_2_sim,
    n_samples   = n_samples_Y,
    n_degs      = n_degs_Y,
    log2FC      = log2FC_Y,
    mean_method = mean_method,
    disp_method = disp_method,
    seed        = (seed + 1000)
  )


  ## --- Calculating interaction truth ---
  fX <- sim_X$feat
  fY <- sim_Y$feat

  true_log2FC <- fY$true_log2FC - fX$true_log2FC
  is_DE <- as.numeric(abs(true_log2FC) > 0)
  sim_features <- data.frame(is_DE = is_DE, true_log2FC = true_log2FC)
  rownames(sim_features) <- rownames(fX) # rownames are features (matrix)


  ## --- Verbose message ---
  if (verbose) {
    fmt_counts <- function(meta, label) {
      tab <- table(meta$condition, dnn = NULL)
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = " ")
    }

    message("\n4-group simulation summary:")
    message(sprintf("%s (X):    N: %-4d  n_DEGs: %-4d  log2FC: %-4.1f  %s", sim_ancestry_X, nrow(sim_X$counts), n_degs_X, log2FC_X, fmt_counts(sim_X$meta, ancestry_X)))
    message(sprintf("%s (Y):    N: %-4d  n_DEGs: %-4d  log2FC: %-4.1f  %s", sim_ancestry_Y, nrow(sim_Y$counts), n_degs_Y, log2FC_Y, fmt_counts(sim_Y$meta, ancestry_Y)))
  }


  ## --- Return ---
  return(
    list(
      X = list(
        matr   = sim_X$matr, 
        meta   = sim_X$meta, 
        feat   = sim_X$feat, 
        plot   = sim_X$in_out_plots
      ),
      Y = list(
        matr   = sim_Y$matr, 
        meta   = sim_Y$meta, 
        feat   = sim_Y$feat, 
        plot   = sim_Y$in_out_plots
      ),
      feat = sim_features
    )
  )
}
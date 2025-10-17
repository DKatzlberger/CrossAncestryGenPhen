#' Simulate ancestry-imbalanced cohorts
#'
#' Draw subsamples from two cohorts (X = majority, Y = minority) to mimic
#' imbalanced ancestry designs. Within each ancestry, samples are split
#' into two groups (`g_col`) at a specified ratio.
#'
#' @param X Numeric matrix or data.frame of features for cohort X (majority ancestry).
#' @param Y Numeric matrix or data.frame of features for cohort Y (minority ancestry).
#' @param MX Data.frame of metadata for X; must align with rows of X.
#' @param MY Data.frame of metadata for Y; must align with rows of Y.
#' @param g_col Column in metadata giving the two-level group split (e.g. case/control).
#' @param a_col Column in metadata giving the ancestry label (one unique value in each of MX, MY).
#' @param total_samples Total number of samples drawn across both ancestries.
#' @param between_ratio Ratio of majority:minority ancestry (e.g. 2 = 2:1).
#' @param within_major_ratio Ratio of group1:group2 within majority ancestry.
#' @param within_minor_ratio Ratio of group1:group2 within minority ancestry.
#' @param seed Optional random seed.
#' @param replace Logical; whether to sample with replacement.
#' @param verbose Logical; print summary messages.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{X}{List with `matr` (subsampled features) and `meta` (aligned metadata) for majority ancestry.}
#'   \item{Y}{List with `matr` (subsampled features) and `meta` (aligned metadata) for minority ancestry.}
#' }
#'
#' @export
sim_imbalanced_ancestry <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col,
  total_samples,
  between_ratio,
  within_major_ratio,
  within_minor_ratio,
  seed = NULL,
  replace = FALSE,
  verbose = TRUE
) {
  ## --- Input data structure check ---
  assert_input(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY, 
    g_col = g_col, 
    a_col = a_col,
    .fun = "sim_imbalanced_ancestry"
  )

  ## --- Seed ---
  if (!is.null(seed)) set.seed(seed)

  ## --- Global sample check ---
  if (!replace) {
    available_total <- nrow(MX) + nrow(MY)
    if (available_total < total_samples) {
      stop(sprintf(
        "[sim_imbalance_ancestry] Not enough samples overall: requested %d, available %d",
        total_samples, available_total
      ))
    }
  }


  ## --- Ancestry levels ---
  a_1 <- unique(MX[[a_col]])
  a_2 <- unique(MY[[a_col]])


  ## --- Helper for sampling within ancestry ---
  sample_two_way <- function(meta, matr, g_col, n, ratio, replace = FALSE) {
    levs <- levels(meta[[g_col]])
    stopifnot(length(levs) == 2)

    r <- as.numeric(ratio)
    p <- r / (r + 1)
    n1_target <- floor(n * p)
    n2_target <- n - n1_target

    i1_all <- which(meta[[g_col]] == levs[1])
    i2_all <- which(meta[[g_col]] == levs[2])
    ## strict check if replace = FALSE
    if (!replace) {
      if (length(i1_all) < n1_target) {
        stop(sprintf("[sim_imbalance_ancestry] Not enough samples in group %s: requested %d, available %d",
                    levs[1], n1_target, length(i1_all)))
      }
      if (length(i2_all) < n2_target) {
        stop(sprintf("[sim_imbalance_ancestry] Not enough samples in group %s: requested %d, available %d",
                    levs[2], n2_target, length(i2_all)))
      }
    }

    n1 <- if (replace) n1_target else min(n1_target, length(i1_all))
    n2 <- if (replace) n2_target else min(n2_target, length(i2_all))

    i1 <- sample(i1_all, n1, replace = replace)
    i2 <- sample(i2_all, n2, replace = replace)
    idx <- c(i1, i2)

    list(
      meta = meta[idx, , drop = FALSE],
      matr = matr[idx, , drop = FALSE]
    )
  }

  ## --- Split samples ---
  n_major <- round(total_samples * between_ratio / (between_ratio + 1))
  n_minor <- total_samples - n_major

  X_out <- sample_two_way(MX, X, g_col, n_major, within_major_ratio, replace)
  Y_out <- sample_two_way(MY, Y, g_col, n_minor, within_minor_ratio, replace)

  ## --- Verbose message ---
  if (verbose) {

    fmt_counts <- function(M, g_col) {
      if (is.null(M) || !nrow(M)) return("")
      tab <- table(M[[g_col]], dnn = NULL)
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = "  ")
    }

    # Define group labels for clarity
    grp_X <- interaction(X_out$meta[[g_col]], X_out$meta[[a_col]], drop = TRUE)
    grp_Y <- interaction(Y_out$meta[[g_col]], Y_out$meta[[a_col]], drop = TRUE)
    grp   <- c(grp_X, grp_Y)

    message("\nImbalanced ancestry summary:")
    message(sprintf("%-13s %s", "Groups:", paste(unique(grp), collapse = "  ")))

    message(sprintf(
      "Ancestry (X): %-10s N: %-5d  %-18s  features: %-5d",
      a_1,
      nrow(X_out$meta),
      fmt_counts(X_out$meta, g_col),
      ncol(X_out$matr)
    ))

    message(sprintf(
      "Ancestry (Y): %-10s N: %-5d  %-18s  features: %-5d",
      a_2,
      nrow(Y_out$meta),
      fmt_counts(Y_out$meta, g_col),
      ncol(Y_out$matr)
    ))
  }

  # if (verbose) {
  #   fmt_counts <- function(M, g_col) {
  #     tab <- table(M[[g_col]], dnn = NULL)
  #     paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = " ")
  #   }

  #   # Define group labels for clarity
  #   grp_X <- interaction(X_out$meta[[g_col]], X_out$meta[[a_col]], drop = TRUE)
  #   grp_Y <- interaction(Y_out$meta[[g_col]], Y_out$meta[[a_col]], drop = TRUE)
  #   grp   <- c(grp_X, grp_Y)

  #   message("\nImbalanced ancestry summary:")
  #   message(sprintf("Groups:    %s", paste(unique(grp), collapse = "  ")))
  #   message(sprintf("%s (X):    N = %-4d %s features: %-4d", a_1, nrow(X_out$meta), fmt_counts(X_out$meta, g_col), ncol(X_out$matr)))
  #   message(sprintf("%s (Y):    N = %-4d %s features: %-4d", a_2, nrow(Y_out$meta), fmt_counts(Y_out$meta, g_col), ncol(Y_out$matr)))
  # }


  ## --- Return ---
  list(
    X = list(
      matr = X_out$matr,
      meta = X_out$meta
    ),
    Y = list(
      matr = Y_out$matr, 
      meta = Y_out$meta
    )
  )
}

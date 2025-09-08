#' Simulate imbalanced ancestry sampling across two cohorts
#'
#' Subsample two data matrices (X, Y) and their metadata (MX, MY) to
#' mimic an ancestry-imbalanced cohort with controllable between- and
#' within-ancestry group splits.
#'
#' X and MX must have the same number of rows; Y and MY must match as
#' well. The ancestry of each block is read from `a_col` in MX and MY;
#' each must contain exactly one non-NA, single label, and the two
#' labels must differ.
#'
#' The within-ancestry split uses `g_col`, which must have exactly two
#' non-NA labels in each block. Labels are inferred alphabetically, and
#' a ratio `r` means "`level1:level2 = r:1`".
#'
#' Targets are snapped so that the within-ancestry ratio is exact
#' whenever possible (using floor for the first group, leftover to the
#' second). If `replace = FALSE` and there are not enough rows, the
#' function draws what is available. A single compact summary of the
#' resulting subsets is printed when `verbose = TRUE`.
#'
#' @param X Numeric matrix or data.frame of features for cohort X; rows are samples and must align with `MX`.
#' @param Y Numeric matrix or data.frame of features for cohort Y; rows are samples and must align with `MY`.
#' @param MX Data.frame with metadata for `X`.
#' @param MY Data.frame with metadata for `Y`.
#' @param g_col Name of the metadata column with the two-level group split used inside each ancestry (e.g., case/control).
#' @param a_col Name of the metadata column holding the ancestry label.
#' @param majority Character, which block is the majority ancestry; one of `c("X","Y")`.
#' @param total_samples Integer, total number of samples drawn across both ancestries.
#' @param between_ratio Positive number, majority:minority ancestry ratio. For example `2` means 2:1 between ancestries.
#' @param within_major_ratio Positive number, ratio for the two `g_col` levels inside the majority ancestry (`level1:level2 = r:1`).
#' @param within_minor_ratio Positive number, ratio for the two `g_col` levels inside the minority ancestry (`level1:level2 = r:1`).
#' @param seed Optional integer; if supplied, sets the RNG seed.
#' @param replace Logical; sample with replacement within each block.
#' @param verbose Logical; if `TRUE`, print a final compact summary.
#'
#' @return A list with four elements:
#'   \describe{
#'     \item{X}{Subset of the feature matrix for the first return slot.}
#'     \item{Y}{Subset of the feature matrix for the second return slot.}
#'     \item{MX}{Subset of metadata aligned to `X`.}
#'     \item{MY}{Subset of metadata aligned to `Y`.}
#'   }
#'   Which ancestry lands in `X` vs `Y` depends on `majority`: the
#'   majority ancestry is always returned first.
#'
#' @export
sim_imbalanced_ancestry <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col,
  majority = c("X","Y"),
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
    a_col = a_col
  )

  ## --- Set the seed ---
  if (!is.null(seed)) set.seed(seed)

  ## --- Unique ancestries ---
  uX <- unique(na.omit(MX[[a_col]]))
  uY <- unique(na.omit(MY[[a_col]]))

  ## --- Match arg ---
  majority <- match.arg(majority)

  # --- Assign majority/minority ---
  if (majority == "X") {
    M_major <- MX; D_major <- X; maj_name <- as.character(uX)
    M_minor <- MY; D_minor <- Y; min_name <- as.character(uY)
    x_is_major <- TRUE
  } else {
    M_major <- MY; D_major <- Y; maj_name <- as.character(uY)
    M_minor <- MX; D_minor <- X; min_name <- as.character(uX)
    x_is_major <- FALSE
  }

  # helper: sample n indices from two groups in g_col at ratio (g1:g2)
  sample_two_way <- function(meta, data, g_col, n, ratio,
                             replace = FALSE) {
    g_vals <- as.character(meta[[g_col]])
    levs <- sort(unique(na.omit(g_vals)))
    if (length(levs) != 2) {
      stop("g_col '", g_col, "' must have exactly 2 non-NA labels; found: ",
           paste(levs, collapse = ", "))
    }

    r <- as.numeric(ratio)
    if (length(r) != 1 || !is.finite(r) || r <= 0)
      stop("ratio must be a positive number")

    # easy fix: floor for group1, leftover to group2
    p <- r / (r + 1)
    n1_target <- floor(n * p)
    n2_target <- n - n1_target

    i1_all <- which(g_vals == levs[1])
    i2_all <- which(g_vals == levs[2])

    n1 <- if (replace) n1_target else min(n1_target, length(i1_all))
    n2 <- if (replace) n2_target else min(n2_target, length(i2_all))

    i1 <- if (length(i1_all)) sample(i1_all, n1, replace = replace) else integer(0)
    i2 <- if (length(i2_all)) sample(i2_all, n2, replace = replace) else integer(0)
    idx <- c(i1, i2)

    def1 <- max(0L, n1_target - n1)
    def2 <- max(0L, n2_target - n2)

    if (!replace && (def1 > 0L || def2 > 0L)) {
      stop("Insufficient samples for requested split in '", g_col, "'.",
              call. = FALSE)
    }

    list(
      meta = meta[idx, , drop = FALSE],
      data = data[idx, , drop = FALSE]
    )
  }

  # split total samples by ancestry (majority:minority = between_ratio:1)
  n_major <- round(total_samples * between_ratio / (between_ratio + 1))
  n_minor <- total_samples - n_major

  # sample within each ancestry with its own internal ratio
  maj  <- sample_two_way(M_major, D_major, g_col, n_major, within_major_ratio, replace)
  minr <- sample_two_way(M_minor, D_minor, g_col, n_minor, within_minor_ratio, replace)

  ## --- Final compact verbose summary ---
  if (verbose) {
    fmt_counts <- function(M, g_col) {
      if (is.null(M) || !nrow(M)) return("")
      tab <- table(M[[g_col]], dnn = NULL)
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = " ")
    }

    if (x_is_major) {
      X_name <- maj_name; X_meta <- maj$meta
      Y_name <- min_name; Y_meta <- minr$meta
    } else {
      X_name <- min_name; X_meta <- minr$meta
      Y_name <- maj_name; Y_meta <- maj$meta
    }

    message("\nImbalanced ancestry summary:")
    message(sprintf("%s (X):    N: %-4d %s", X_name, nrow(X_meta), fmt_counts(X_meta, g_col)))
    message(sprintf("%s (Y):    N: %-4d %s", Y_name, nrow(Y_meta), fmt_counts(Y_meta, g_col)))
  }

  ## --- Return ---
  if (x_is_major) {
    list(
      X = list(counts = maj$data,  meta = maj$meta),
      Y = list(counts = minr$data, meta = minr$meta)
    )
  } else {
    list(
      X = list(counts = minr$data, meta = minr$meta),
      Y = list(counts = maj$data,  meta = maj$meta)
    )
  }
}

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
#' Targets are rounded. If `replace = FALSE` and there are not enough
#' rows, the function draws what is available. Per-block split
#' messages are printed only when `verbose = TRUE`. A warning is
#' emitted on shortfall when sampling without replacement.
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
#' @param verbose Logical; if `TRUE`, print per-block split messages.
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
#' @examples
#' set.seed(1)
#' X  <- matrix(rnorm(1000), nrow = 100)
#' Y  <- matrix(rnorm( 800), nrow =  80)
#' MX <- data.frame(anc = "EUR",
#'                  grp = sample(c("A","B"), 100, TRUE))
#' MY <- data.frame(anc = "AFR",
#'                  grp = sample(c("A","B"),  80, TRUE))
#' out <- sim_imbalanced_ancestry(
#'   X, Y, MX, MY,
#'   g_col = "grp",
#'   a_col = "anc",
#'   majority = "X",
#'   total_samples = 120,
#'   between_ratio = 2,
#'   within_major_ratio = 3,
#'   within_minor_ratio = 1,
#'   seed = 123,
#'   verbose = TRUE
#' )
#' lapply(out, nrow)
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
  total_samples = 100,
  between_ratio = 1,
  within_major_ratio = 1,
  within_minor_ratio = 1,
  seed = NULL,
  replace = FALSE,
  verbose = FALSE
) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(nrow(X) == nrow(MX), nrow(Y) == nrow(MY))
  majority <- match.arg(majority)

  # column checks
  if (!(g_col %in% names(MX)) || !(g_col %in% names(MY)))
    stop("g_col must exist in both MX and MY.")
  if (!(a_col %in% names(MX)) || !(a_col %in% names(MY)))
    stop("a_col must exist in both MX and MY.")

  # ancestry labels per block (must be single, distinct values)
  uX <- unique(na.omit(MX[[a_col]]))
  uY <- unique(na.omit(MY[[a_col]]))
  if (length(uX) != 1L)
    stop("In MX, a_col must have exactly one non-NA value.")
  if (length(uY) != 1L)
    stop("In MY, a_col must have exactly one non-NA value.")
  if (identical(as.character(uX), as.character(uY)))
    stop("a_col must differ between MX and MY (it is the split column).")

  # assign majority/minority metas + names
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
                             block_name = "", replace = FALSE,
                             verbose = TRUE) {
    g_vals <- as.character(meta[[g_col]])
    # infer labels from g_col (alphabetical)
    levs <- sort(unique(na.omit(g_vals)))
    if (length(levs) != 2) {
      stop("In block ", block_name, ", g_col '", g_col,
           "' must have exactly 2 non-NA labels; found: ",
           paste(levs, collapse = ", "))
    }

    # targets before availability check
    r <- as.numeric(ratio)
    if (length(r) != 1 || !is.finite(r) || r <= 0)
      stop("ratio must be a positive number")
    n1_target <- round(n * r / (r + 1))
    n2_target <- n - n1_target

    # available indices (compare on character labels)
    i1_all <- which(g_vals == levs[1])
    i2_all <- which(g_vals == levs[2])

    # actual draws (respect availability if replace = FALSE)
    n1 <- if (replace) n1_target else min(n1_target, length(i1_all))
    n2 <- if (replace) n2_target else min(n2_target, length(i2_all))

    i1 <- if (length(i1_all)) sample(i1_all, n1, replace = replace) else integer(0)
    i2 <- if (length(i2_all)) sample(i2_all, n2, replace = replace) else integer(0)
    idx <- c(i1, i2)

    # deficits
    def1 <- max(0L, n1_target - n1)
    def2 <- max(0L, n2_target - n2)

    # print a message only when verbose
    if (isTRUE(verbose)) {
      message(paste0(
        "[", block_name, "] Split in ", g_col, ": ",
        "Target {", levs[1], ":", n1_target, "; ",
        levs[2], ":", n2_target, "}. ",
        "Drew {", levs[1], ":", n1, "; ",
        levs[2], ":", n2, "}. ",
        "Missing {", levs[1], ":", def1, "; ",
        levs[2], ":", def2, "}."
      ))
    }

    # warn only if there's a shortfall without replacement
    if (!replace && (def1 > 0L || def2 > 0L)) {
      warning(paste0("[", block_name,
                     "] Insufficient samples for requested split."),
              call. = FALSE)
    }

    # return the subset matrices/metadata
    list(
      meta = meta[idx, , drop = FALSE],
      data = data[idx, , drop = FALSE]
    )
  }

  # split total samples by ancestry (majority:minority = between_ratio:1)
  n_major <- round(total_samples * between_ratio / (between_ratio + 1))
  n_minor <- total_samples - n_major

  # sample within each ancestry with its own internal ratio
  maj <- sample_two_way(
    meta = M_major, data = D_major, g_col = g_col, n = n_major,
    ratio = within_major_ratio,
    block_name = paste0("Majority ", maj_name),
    replace = replace, verbose = verbose
  )
  minr <- sample_two_way(
    meta = M_minor, data = D_minor, g_col = g_col, n = n_minor,
    ratio = within_minor_ratio,
    block_name = paste0("Minority ", min_name),
    replace = replace, verbose = verbose
  )

  # return ONLY expression/count matrices (X,Y) and their metadata (MX,MY)
  if (x_is_major) {
    return(list(
      X = maj$data, # matrix with rownames = samples
      Y = minr$data,
      MX = maj$meta,
      MY = minr$meta
    ))
  } else {
    return(list(
      X = minr$data,
      Y = maj$data,
      MX = minr$meta,
      MY = maj$meta
    ))
  }
}

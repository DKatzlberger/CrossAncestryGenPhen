#' Split Expression and Metadata into Reference (R), Subset (X), and Inference (Y) Sets 
#'
#' Performs stratified sampling of an overrepresented group (X, e.g. EUR) to match the distribution of 
#' an underrepresented group (Y, e.g. AFR) based on a grouping variable (e.g., condition).
#'
#' @param X Numeric matrix or data.frame of features for cohort X; rows are samples and must align with MX.
#' @param Y Numeric matrix or data.frame of features for cohort Y; rows are samples and must align with MY.
#' @param MX Data.frame with metadata for X.
#' @param MY Data.frame with metadata for Y.
#' @param g_col Name of the metadata column holding the stratification label.
#' @param a_col Name of the metadata column holding the ancestry label.
#' @param match Logical, whether to subset both X and Y to have 'a x g' balance.
#' @param seed Optional numeric seed for reproducibility of sampling.
#' @param verbose Logical, whether to print messages.
#'
#' @return A list with the following elements (all matrices with rownames):
#'   \item{RX}{Reference set: remaining X after subsampling}
#'   \item{RY}{Reference set: remaining Y after subsampling (if match = TRUE)}
#'   \item{X}{Subset set: subsampled X matching Y}
#'   \item{Y}{Inference set: subsampled Y (if match = TRUE) or full Y}
#'   \item{strata_info}{list with usable/missing/insufficient strata}
#' 
#' @export
split_stratified_ancestry_sets <- function(
  X, 
  Y, 
  MX, 
  MY,
  g_col, 
  a_col,
  match = FALSE,
  seed = NULL,
  verbose = TRUE
) {
    
  ## --- Seed ---
  if (!is.null(seed)) set.seed(seed)

  ## --- Input checks ---
  assert_input(
    X = X, 
    Y = Y,
    MX = MX, 
    MY = MY,
    g_col = g_col, 
    a_col = a_col
  )

  ## --- Factor setup ---
  a_1 <- unique(MX[[a_col]])
  a_2 <- unique(MY[[a_col]])

  g_levels <- levels(MX[[g_col]])
  if (length(g_levels) != 2 || length(unique(c(a_1, a_2))) != 2) {
    stop("[split_stratified_ancestry_sets] Function supports only 2x2 designs (two levels in g_col x two levels a_col).")
  }

  # KORREKTUR: Erzwinge exakt dieselben Faktorstufen für beide Kohorten
  vec_g_X <- factor(MX[[g_col]], levels = g_levels)
  vec_g_Y <- factor(MY[[g_col]], levels = g_levels)

  ## --- Target Count Calculation ---
  count_X <- table(vec_g_X)
  count_Y <- table(vec_g_Y)

  if (match) {
    # Berechne das absolute Minimum über alle Gruppen hinweg
    min_overall   <- min(c(count_X, count_Y))
    target_counts <- setNames(rep(min_overall, length(g_levels)), g_levels)
  } else {
    target_counts <- count_Y
  }

  strata_names <- names(target_counts)

  ## --- Feasibility check ---
  insufficient <- strata_names[target_counts[strata_names] > count_X[strata_names]]
  missing      <- setdiff(strata_names, names(count_X))

  if (length(missing) > 0 || length(insufficient) > 0) {
    stop("[split_stratified_ancestry_sets] X cannot fulfill the requested distribution layout.\n",
         "Missing strata: ", paste(missing, collapse = ", "), "\n",
         "Insufficient strata: ", paste(insufficient, collapse = ", "))
  }

  ## --- Process Y (Inference & Remaining RY) ---
  ids_Y <- rownames(Y)
  
  if (match) {
    sampled_ids_Y <- vector("list", length(strata_names))
    for (i in seq_along(strata_names)) {
      stratum <- strata_names[i]
      idx <- which(vec_g_Y == stratum)
      # target_counts[stratum] zieht jetzt garantiert die exakt korrekte Anzahl für diesen Namen
      sampled_ids_Y[[i]] <- ids_Y[sample(idx, size = target_counts[stratum], replace = FALSE)]
    }
    sampled_ids_Y <- unlist(sampled_ids_Y, use.names = FALSE)
    
    mask_Y_subset <- ids_Y %in% sampled_ids_Y
    
    Y_matr  <- Y[mask_Y_subset, , drop = FALSE]
    Y_meta  <- MY[mask_Y_subset, , drop = FALSE]
    RY_matr <- Y[!mask_Y_subset, , drop = FALSE]
    RY_meta <- MY[!mask_Y_subset, , drop = FALSE]
  } else {
    Y_matr  <- Y
    Y_meta  <- MY
    RY_matr <- NULL
    RY_meta <- NULL
  }

  ## --- Process X (Subset X & Remaining RX) ---
  ids_X <- rownames(X)
  sampled_ids_X <- vector("list", length(strata_names))
  for (i in seq_along(strata_names)) {
    stratum <- strata_names[i]
    idx <- which(vec_g_X == stratum)
    sampled_ids_X[[i]] <- ids_X[sample(idx, size = target_counts[stratum], replace = FALSE)]
  }
  sampled_ids_X <- unlist(sampled_ids_X, use.names = FALSE)

  mask_X_subset <- ids_X %in% sampled_ids_X

  X_matr  <- X[mask_X_subset, , drop = FALSE]
  X_meta  <- MX[mask_X_subset, , drop = FALSE]
  RX_matr <- X[!mask_X_subset, , drop = FALSE]
  RX_meta <- MX[!mask_X_subset, , drop = FALSE]

  ## --- Verbose summary ---
  if (verbose) {
    fmt_counts <- function(M_sub, g_col) {
      if (is.null(M_sub) || nrow(M_sub) == 0) return("N/A")
      # Gewährleistet konsistente Level-Reihenfolge in der Konsolen-Ausgabe
      tab <- table(factor(M_sub[[g_col]], levels = g_levels))
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = " ")
    }

    message("\nStratified split:")
    if (match) {
      message("Enforcing 'a_col x g_col' balance.")
      message(sprintf("%-20s  N: %-4d %s features: %-4d", sprintf("Remaining RX %-8s", paste0("(", a_1, ")")), nrow(RX_matr), fmt_counts(RX_meta, g_col), ncol(RX_matr)))
      message(sprintf("%-20s  N: %-4d %s features: %-4d", sprintf("Remaining RY %-8s", paste0("(", a_2, ")")), nrow(RY_matr), fmt_counts(RY_meta, g_col), ncol(RY_matr)))
      message(sprintf("%-20s  N: %-4d %s features: %-4d", sprintf("Subset    SX %-8s", paste0("(", a_1, ")")), nrow(X_matr), fmt_counts(X_meta, g_col), ncol(X_matr)))
      message(sprintf("%-20s  N: %-4d %s features: %-4d", sprintf("Subset    SY %-8s", paste0("(", a_2, ")")), nrow(Y_matr), fmt_counts(Y_meta, g_col), ncol(Y_matr)))
    } else {
      message("Enforcing 'a_col' balance.")
      message(sprintf("%-20s  N: %-4d %s features: %-4d", sprintf("Remaining RX %-8s", paste0("(", a_1, ")")), nrow(RX_matr), fmt_counts(RX_meta, g_col), ncol(RX_matr)))
      message(sprintf("%-20s  N: %-4d %s features: %-4d", sprintf("Subset    SX %-8s", paste0("(", a_1, ")")), nrow(X_matr), fmt_counts(X_meta, g_col), ncol(X_matr)))
      message(sprintf("%-20s  N: %-4d %s features: %-4d", sprintf("Original  Y  %-8s", paste0("(", a_2, ")")), nrow(Y_matr), fmt_counts(Y_meta, g_col), ncol(Y_matr)))
    }
  }

  ## --- Return ---
  return(
    list(
      RX = list(matr = RX_matr, meta = RX_meta, ids = rownames(RX_matr)),
      RY = list(matr = RY_matr, meta = RY_meta, ids = if(!is.null(RY_matr)) rownames(RY_matr) else NULL),
      X  = list(matr = X_matr,  meta = X_meta,  ids = rownames(X_matr)),
      Y  = list(matr = Y_matr,  meta = Y_meta,  ids = rownames(Y_matr)),
      strata_info = list(usable = strata_names, missing = missing, insufficient = insufficient)
    )
  )
}

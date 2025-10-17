#' Subset metadata and expression by ancestry and phenotype
#'
#' Filters metadata (`M`) and matching matrix (`X`) to two given ancestry
#' levels and two phenotype levels, then splits into X/Y groups.
#'
#' @param X Matrix of features (rows = samples, cols = features).
#' @param M Data frame of sample metadata (rownames = sample IDs).
#' @param g_col Column in `M` with phenotype labels.
#' @param a_col Column in `M` with ancestry labels.
#' @param g_levels Character(2), phenotype levels to keep (order sets factor).
#' @param a_levels Character(2), ancestry levels to keep (`X`, `Y`).
#' @param verbose Logical, print summary (default TRUE).
#'
#' @return List with elements `X` and `Y`, each a list of:
#'   \itemize{
#'     \item \code{meta}: filtered metadata
#'     \item \code{counts}: filtered matrix
#'   }
#'
#' @examples
#' \dontrun{
#' subset_phenotype_ancestry(
#'   X = study$mrna, M = study$meta,
#'   g_col = "subtype", a_col = "ancestry",
#'   g_levels = c("Basal", "LumA"),
#'   a_levels = c("EUR", "AFR")
#' )
#' }
#' @export
filter_phenotype_ancestry <- function(
  X,
  M,
  g_col,
  a_col,
  g_levels,
  a_levels,
  covariates = NULL,
  plot = TRUE,
  verbose = TRUE
){

  ## --- Input data structure check ---
  ## Different from how normal checked becuase at the beginning of the pipeline
  if (length(g_levels) != 2 || length(a_levels) != 2) {
    stop("[filter_phenotype_ancestry] Function supports only 2x2 designs (two levels in g_col × two levels a_col).")
  }

  stopifnot(is.matrix(X), !is.null(rownames(X)))
  stopifnot(is.data.frame(M), !is.null(rownames(M)))
  stopifnot(g_col %in% colnames(M), a_col %in% colnames(M))

  if (!identical(rownames(X), rownames(M))) {
    stop("Row names of `M` and `X` must match exactly (same sampleIDs).")
  }


  ## --- Filter meta by columns ---
  meta_cols <- unique(c(g_col, a_col, covariates))
  meta_cols <- meta_cols[meta_cols %in% colnames(M)]

  if (length(meta_cols) == 0) {
    stop("[filter_phenotype_ancestry] No valid metadata columns found (check g_col, a_col, and covariates).")
  }


  ## --- Filter meta rows (M) ---
  keep   <- M[[a_col]] %in% a_levels & M[[g_col]] %in% g_levels
  M_filt <- M[keep, , drop = FALSE]


  ## --- Filter counts (X) ---
  X_filt <- X[rownames(M_filt), , drop = FALSE]


  ## --- Split by ancestry ---
  split_by_ancestry <- function(a){

    # Subset by ancestry
    m_sub <- M_filt[M_filt[[a_col]] == a, , drop = FALSE]

    if (nrow(m_sub) == 0L) {
      stop(sprintf("No samples found for ancestry level '%s'.", a))
    }
    
    # Enforce factor levels for g_col
    if (g_col %in% colnames(m_sub)) {
      m_sub[[g_col]] <- factor(m_sub[[g_col]], levels = g_levels)
    }
    
    ids <- rownames(m_sub)
    m_sub <- m_sub[ids, , drop = FALSE]
    x_sub <- X[ids, , drop = FALSE]
    
    # Run your existing assertion helper
    assert_input(
      X  = x_sub,
      MX = m_sub,
      g_col = g_col,
      a_col = a_col,
      .fun  = "subset_phenotype_ancestry"
    )

    
    list(
      meta = m_sub, 
      matr = x_sub,
      ancestry = a
    )
  }

  ## --- Apply function ---
  X_out <- split_by_ancestry(a_levels[1])
  Y_out <- split_by_ancestry(a_levels[2])


  ## --- NA check for each ancestry ---
  na_X <- colSums(is.na(X_out$meta[, meta_cols, drop = FALSE]))
  na_Y <- colSums(is.na(Y_out$meta[, meta_cols, drop = FALSE]))

  total_na_X <- sum(na_X)
  total_na_Y <- sum(na_Y)

  if (total_na_X > 0 || total_na_Y > 0) {

    fmt_na <- function(na_vec) {
      if (is.null(na_vec) || !length(na_vec)) return("")
      paste(sprintf("%s: %-4d", names(na_vec), as.integer(na_vec)), collapse = "  ")
    }

    msg_X <- fmt_na(na_X)
    msg_Y <- fmt_na(na_Y)

    message("\nNAs detected:")

    message(sprintf(
      "Ancestry (X): %-10s  N: %-5d  %-40s  | Total NAs: %-6d",
      X_out$ancestry,
      nrow(X_out$meta),
      msg_X,
      total_na_X
    ))
    
    message(sprintf(
      "Ancestry (Y): %-10s  N: %-5d  %-40s  | Total NAs: %-6d",
      Y_out$ancestry,
      nrow(Y_out$meta),
      msg_Y,
      total_na_Y
    ))

    stop()
  }


  ## --- Imbalace plot ---
  p <- plot_imbalanced_groups(
    MX = X_out$meta,
    MY = Y_out$meta,
    x_var = a_col,
    fill_var = g_col,
    title = NULL,
    x_label = NULL,
    y_label = NULL
  )

  if (plot){
    print(p)
  }


  ## --- Verbose message ----
  if (verbose) {

    fmt_counts <- function(M, g_col) {
      if (is.null(M) || !nrow(M)) return("")
      tab <- table(M[[g_col]], dnn = NULL)
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = "  ")
    }

    # Defined groups
    grp_X <- interaction(X_out$meta[[g_col]], X_out$meta[[a_col]], drop = TRUE)
    grp_Y <- interaction(Y_out$meta[[g_col]], Y_out$meta[[a_col]], drop = TRUE)
    grp   <- c(grp_X, grp_Y)

    message("\nSubset phenotype ancestry summary:")
    message(sprintf("%-13s %s", "Groups:", paste(unique(grp), collapse = "  ")))

    message(sprintf(
      "Ancestry (X): %-10s N: %-5d  %-18s  features: %-5d",
      X_out$ancestry,
      nrow(X_out$meta),
      fmt_counts(X_out$meta, g_col),
      ncol(X_out$matr)
    ))

    message(sprintf(
      "Ancestry (Y): %-10s N: %-5d  %-18s  features: %-5d",
      Y_out$ancestry,
      nrow(Y_out$meta),
      fmt_counts(Y_out$meta, g_col),
      ncol(Y_out$matr)
    ))
  }

  # if (verbose) {
  #   fmt_counts <- function(M, g_col) {
  #     if (is.null(M) || !nrow(M)) return("")
  #     tab <- table(M[[g_col]], dnn = NULL)
  #     paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = " ")
  #   }


  #   # Defined groups
  #   grp_X <- interaction(X_out$meta[[g_col]], X_out$meta[[a_col]], drop = TRUE)
  #   grp_Y <- interaction(Y_out$meta[[g_col]], Y_out$meta[[a_col]], drop = TRUE)
  #   grp   <- c(grp_X, grp_Y)

  #   message("\nSubset phenotype ancestry summary:")
  #   message(sprintf("Groups:     %s", paste(unique(grp), collapse = "  ")))
  #   message(sprintf("%s (X):    N: %-4d %s features: %-4d", X_out$ancestry, nrow(X_out$meta), fmt_counts(X_out$meta, g_col), ncol(X_out$matr)))
  #   message(sprintf("%s (Y):    N: %-4d %s features: %-4d", Y_out$ancestry, nrow(Y_out$meta), fmt_counts(Y_out$meta, g_col), ncol(Y_out$matr)))
  # }


  ## --- Return ---
  out <- list(
    X = X_out,
    Y = Y_out
  )

  # Add the plot
  out$plot <- p

  return(out)
}
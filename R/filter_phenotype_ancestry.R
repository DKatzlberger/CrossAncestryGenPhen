#' Subset metadata and expression by ancestry and phenotype
#'
#' Filters metadata (`M`) and a matching feature matrix (`X`) to two selected
#' ancestry levels and two phenotype levels, then splits the data into two
#' groups (X/Y) based on ancestry. Optionally removes samples with missing
#' values in required metadata columns.
#'
#' @param X A numeric matrix of features (rows = samples, columns = features).
#' @param M A data frame containing sample metadata (must have rownames matching `rownames(X)`).
#' @param g_col A character string specifying the metadata column containing phenotype labels.
#' @param a_col A character string specifying the metadata column containing ancestry labels.
#' @param g_levels A character vector of length 2 giving the phenotype levels to retain (order determines reference vs. comparison).
#' @param a_levels A character vector of length 2 giving the ancestry levels to retain (first becomes group X, second group Y).
#' @param covariates Optional character vector of additional metadata columns to include during NA checking (e.g., age).
#' @param omit_na Logical; if `TRUE`, missing values in any of the required metadata columns (`g_col`, `a_col`, `covariates`) are removed.
#' @param auto_cast Logical; if `TRUE`, covariates are automatically cast into classes.
#' @param plot Logical; if `TRUE` (default), prints the imbalance plot showing phenotype frequencies across ancestry groups.
#' @param verbose Logical; if `TRUE` (default), prints detailed progress and summary information.
#'
#' @export
filter_phenotype_ancestry <- function(
  X,
  M,
  g_col,
  a_col,
  g_levels,
  a_levels,
  covariates = NULL,
  omit_na = FALSE,
  auto_cast = FALSE,
  plot = TRUE,
  verbose = TRUE
){

  ## --- Helper: Split by ancestry ---
  split_by_ancestry <- function(a, M_current, X_current, g_col, a_col, g_levels) {

    # Subset metadata
    m_sub <- M_current[M_current[[a_col]] == a, , drop = FALSE]
    if (nrow(m_sub) == 0L) stop(sprintf("No samples found for ancestry level '%s'.", a))

    # Enforce phenotype factor order
    m_sub[[g_col]] <- factor(
      as.character(m_sub[[g_col]]),
      levels = g_levels
    )

    # Subset X accordingly
    ids <- rownames(m_sub)
    x_sub <- X_current[ids, , drop = FALSE]

    # Validate
    assert_input(
      X  = x_sub,
      MX = m_sub,
      g_col = g_col,
      a_col = a_col,
      .fun  = "subset_phenotype_ancestry"
    )

    # Return clean structured output
    list(
      meta     = m_sub,
      matr     = x_sub,
      ancestry = a
    )
  }


  ## --- Helper: Covariate detection ---
  cast_covariate_columns <- function(M, cols, continuous_threshold = 3) {

    # regex for integer, decimal, or scientific notation
    continuous_regex <- "^\\s*[+-]?((\\d+\\.?\\d*)|(\\d*\\.\\d+))(e[+-]?\\d+)?\\s*$"

    for (cc in cols) {

      v <- M[[cc]]

      # 1. Never modify factors
      if (is.factor(v)) next
      # 2. Never modify numeric or logical
      if (is.numeric(v) || is.logical(v)) next
      # 3. Only characters can be continuous-like
      if (is.character(v)) {
        numeric_form <- all(grepl(continuous_regex, v[!is.na(v)], ignore.case = TRUE))
        if (numeric_form) {
          v_num <- as.numeric(v)
          nuniq <- length(unique(v_num[!is.na(v_num)]))

          # Only cast if it is truly continuous
          if (nuniq > continuous_threshold) {
            M[[cc]] <- v_num
          }
        }

        next
      }
    }

    return(M)
  }



  ## --- Input data structure check (Different from how normal checked becuase at the beginning of the pipeline) ---
  if (length(g_levels) != 2 || length(a_levels) != 2) {
    stop("[filter_phenotype_ancestry] Function supports only 2x2 designs (two levels in g_col × two levels a_col).")
  }

  stopifnot(is.matrix(X), !is.null(rownames(X)))
  stopifnot(is.data.frame(M), !is.null(rownames(M)))
  stopifnot(g_col %in% colnames(M), a_col %in% colnames(M))
  if (!identical(rownames(X), rownames(M))) stop("Row names of `M` and `X` must match exactly (same sampleIDs).")
  


  ## --- Filter meta by columns ---
  meta_cols <- unique(c(g_col, a_col, covariates))
  meta_cols <- meta_cols[meta_cols %in% colnames(M)]


  if (length(meta_cols) == 0) {
    stop("[filter_phenotype_ancestry] No valid metadata columns found (check g_col, a_col, and covariates).")
  }

  ## --- Validate ancestry and phenotype columns ---
  if (!a_col %in% colnames(M)) {
    stop(sprintf("[filter_phenotype_ancestry] a_col = '%s' is NOT a column in metadata.", a_col))
  }

  if (!g_col %in% colnames(M)) {
    stop(sprintf("[filter_phenotype_ancestry] g_col = '%s' is NOT a column in metadata.", g_col))
  }

  ## --- Validate requested ancestry levels ---
  missing_a <- setdiff(a_levels, unique(M[[a_col]]))
  if (length(missing_a) > 0) {
    stop(sprintf("[filter_phenotype_ancestry] Requested ancestry level(s) NOT found: %s", paste(missing_a, collapse = ", ")))
  }

  ## --- Validate requested phenotype levels ---
  missing_g <- setdiff(g_levels, unique(M[[g_col]]))
  if (length(missing_g) > 0) {
    stop(sprintf("[filter_phenotype_ancestry] Requested phenotype level(s) NOT found: %s", paste(missing_g, collapse = ", ")))
  }


  ## --- Filter meta rows (M) ---
  keep   <- M[[a_col]] %in% a_levels & M[[g_col]] %in% g_levels
  M_filt <- M[keep, , drop = FALSE]


  ## --- Apply split function ---
  X_out <- split_by_ancestry(a_levels[1], M_filt, X, g_col, a_col, g_levels)
  Y_out <- split_by_ancestry(a_levels[2], M_filt, X, g_col, a_col, g_levels)


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

    if (verbose){
      message("\nNAs detected:")

      message(sprintf(
        "Ancestry (X): %-10s N: %-5d  %s  | Total NAs: %-6d",
        X_out$ancestry,
        nrow(X_out$meta),
        msg_X,
        total_na_X
      ))

      message(sprintf(
        "Ancestry (Y): %-10s N: %-5d  %s  | Total NAs: %-6d",
        Y_out$ancestry,
        nrow(Y_out$meta),
        msg_Y,
        total_na_Y
      ))
    }

    ## --- Drop NA rows if requested ---
    if (omit_na) {
      if (verbose) message("\nomit_na = TRUE → Removing rows with NA in required metadata columns...")

      keep_complete <- complete.cases(M_filt[, meta_cols, drop = FALSE])
      M_filt <- M_filt[keep_complete, , drop = FALSE]

      # Re-factor phenotype
      M_filt[[g_col]] <- factor(M_filt[[g_col]], levels = g_levels)

    } else {
      stop("NAs detected. Use omit_na = TRUE to drop NA rows.")
    }
  }

  ## --- ALWAYS final split, NA or not ---
  X_out <- split_by_ancestry(a_levels[1], M_filt, X, g_col, a_col, g_levels)
  Y_out <- split_by_ancestry(a_levels[2], M_filt, X, g_col, a_col, g_levels)


  ## --- Covariates detected ---
  if (length(meta_cols) > 0) {

    detect_covariate_info <- function(M, cols, continuous_threshold = 3) {

      continuous_regex <- "^\\s*[+-]?((\\d+\\.?\\d*)|(\\d*\\.\\d+))(e[+-]?\\d+)?\\s*$"
      out <- lapply(cols, function(cc) {

        v   <- M[[cc]]
        cls <- class(v)[1]

        # Count unique values
        nuniq <- length(unique(v[!is.na(v)]))

        # Predict casting (same logic as above)
        will_cast <- FALSE

        if (is.character(v)) {

          numeric_form <- all(grepl(continuous_regex, v[!is.na(v)], ignore.case = TRUE))

          if (numeric_form) {
            v_num <- as.numeric(v)
            nuniq_num <- length(unique(v_num[!is.na(v_num)]))
            will_cast <- nuniq_num > continuous_threshold
          }

        } else if (is.numeric(v)) {

          # numeric columns remain numeric, not "cast"
          will_cast <- FALSE

        } else {
          # factors, logicals, etc. never cast
          will_cast <- FALSE
        }

        list(
          class     = cls,
          unique_n  = nuniq,
          cast_to   = if (will_cast) "continuous" else cls
        )
      })

      names(out) <- cols
      return(out)
    }

    info_X <- detect_covariate_info(X_out$meta, meta_cols)
    info_Y <- detect_covariate_info(Y_out$meta, meta_cols)

    if (verbose) {
      message("\nCovariates detected:")

      fmt_covline <- function(info) {
        paste(
          sapply(names(info), function(cc) {
            e <- info[[cc]]
            sprintf(
              "%s: %s -> %s",
              cc,
              e$class,
              e$cast_to
            )
          }),
          collapse = "  "
        )
      }

      message(sprintf(
        "Ancestry (X): %-10s N: %-5d  %s",
        X_out$ancestry, nrow(X_out$meta), fmt_covline(info_X)
      ))

      message(sprintf(
        "Ancestry (Y): %-10s N: %-5d  %s",
        Y_out$ancestry, nrow(Y_out$meta), fmt_covline(info_Y)
      ))
    }
  }


  ## --- Auto cast class for covariates ---
  if (auto_cast) {
    if (verbose) message("\nauto_cast = TRUE → Converting class in required metadata columns...")

    # X group
    cast_X <- if (length(meta_cols) > 0) {
      cast_covariate_columns(X_out$meta, meta_cols)
    } else X_out$meta
    X_out$meta <- cast_X

    # Y group
    cast_Y <- if (length(meta_cols) > 0) {
      cast_covariate_columns(Y_out$meta, meta_cols)
    } else Y_out$meta
    Y_out$meta <- cast_Y
  }


  ## --- Check that covariates match (only if covariates exist) ---
  cov_info <- NULL
  if (length(meta_cols) > 0) {
    cov_info <- sapply(meta_cols, function(cc) {
      class_X <- class(X_out$meta[[cc]])[1]
      class_Y <- class(Y_out$meta[[cc]])[1]

      if (!identical(class_X, class_Y)) {
        stop(sprintf(
          "Covariate '%s' has mismatched classes between ancestry groups: X = %s, Y = %s",
          cc, class_X, class_Y
        ))
      }
      class_X
    }, USE.NAMES = TRUE)
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

    
    # Define covariates
    if (length(meta_cols) > 0) {
      cov_line <- paste(sprintf("%s (%s)", names(cov_info), cov_info), collapse = "  ")
    }

    message("\nSubset phenotype ancestry summary:")
    message(sprintf("%-13s %s", "Groups:", paste(unique(grp), collapse = "  ")))
    if (length(meta_cols) > 0) message(sprintf("%-13s %s", "Covariates:", cov_line))
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
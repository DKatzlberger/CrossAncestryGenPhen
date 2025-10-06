#' DESeq2-Based Interaction Test (Wald)
#'
#' @param X Expression matrix for ancestry X. Rows = samples, columns = genes.
#' @param Y Expression matrix for ancestry Y. Rows = samples, columns = genes.
#' @param MX Metadata for X. Must include group and ancestry columns.
#' @param MY Metadata for Y. Must include group and ancestry columns.
#' @param g_col Name of the column indicating group (factor 1).
#' @param a_col Name of the column indicating ancestry (factor 2).
#' @param covariates Optional vector of covariate column names to adjust for.
#' @param verbose Logical, whether to print summary messages.
#'
#' @return summary stats
#'
#' @export
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results resultsNames counts
#' @importFrom limma makeContrasts
DESeq_interaction_effect <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col,
  covariates = NULL,
  verbose = TRUE
){

  ## --- Input data structure check ---
  assert_input(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY, 
    g_col = g_col, 
    a_col = a_col,
    .fun = "DESeq_interaction_effect"
  )
  

  ## --- Combine expression and metadata ---
  matr <- rbind(X, Y)
  meta <- rbind(MX, MY)

  if (!identical(rownames(matr), rownames(meta))) {
    stop("[DESeq_interaction_effect] Matrix and meta rownames must match exactly.")
  }


  ## --- Factor setup ---
  a_1 <- unique(MX[[a_col]])
  a_2 <- unique(MY[[a_col]])

  meta[[a_col]] <- factor(meta[[a_col]], levels = c(a_1, a_2))

  g_levels <- levels(meta[[g_col]])
  a_levels <- levels(meta[[a_col]])

  if (length(g_levels) != 2 || length(a_levels) != 2) {
    stop("[DESeq_interaction_effect] Function currently supports only 2x2 designs (two levels in g_col × two levels a_col).")
  }


  ## --- Build 4 groups ----
  g_1 <- g_levels[1]; g_2 <- g_levels[2]

  meta[["groups"]] <- factor(
    paste(
        meta[[g_col]], 
        meta[[a_col]], 
        sep = "."
    ),
    levels = c(
      paste(g_1, a_1, sep = "."),
      paste(g_2, a_1, sep = "."),
      paste(g_1, a_2, sep = "."),
      paste(g_2, a_2, sep = ".")
    )
  )


  ## --- Build means model formula (4 groups) ---
  form_str <- paste("~0 + groups")
  if (!is.null(covariates)) {
    form_str <- paste(form_str, "+", paste(covariates, collapse = " + "))
  }
  design <- model.matrix(as.formula(form_str), data = meta)
  colnames(design) <- make.names(colnames(design))


  ## --- Group vs covariates coef ---
  all_coefs   <- colnames(design)
  group_mask  <- grepl("^groups", all_coefs)
  group_coefs <- all_coefs[group_mask]
  covar_coefs <- all_coefs[!group_mask]


  ## --- Clean coef names ---
  group_coefs <- gsub("groups", "", group_coefs)               
  colnames(design)[group_mask] <- group_coefs

  if (!is.null(covariates)) {
    for (cov in covariates) {
      if (is.factor(meta[[cov]]) || is.character(meta[[cov]])) {
        covar_coefs <- sub(paste0("^", cov), "", covar_coefs)
      }
    }
    colnames(design)[!group_mask] <- covar_coefs
  }


  ## --- Define contrasts ---
  cols <- colnames(design)[group_mask]

  # Contrast calculations
  contrast_calculations_pretty <- list(
    baseline_1     = paste(cols[3], "-", cols[1]), # G1.A2 - G1.A1
    baseline_2     = paste(cols[4], "-", cols[2]), # G2.A2 - G2.A1
    relationship_1 = paste(cols[2], "-", cols[1]), # G2.A1 - G1.A1
    relationship_2 = paste(cols[4], "-", cols[3]), # G2.A2 - G1.A2
    interaction    = paste0("(", cols[4], " - ", cols[3], ") - (", cols[2], " - ", cols[1], ")")
  )


  ## --- Verbose message ---
  if (verbose) {
    message("\nGLM summary (DESeq):")
    message(sprintf("Formula:         %s", form_str))
    message(sprintf("Groups:          %s", paste(cols, collapse = "  ")))
    message(sprintf("Baseline:        %s", paste(contrast_calculations_pretty[1:2], collapse = "  ")))
    message(sprintf("Relationship:    %s", paste(contrast_calculations_pretty[3:4], collapse = "  ")))
    message(sprintf("Interaction:     %s", paste(contrast_calculations_pretty[[5]], collapse = "  ")))
  }


  ## --- DESeq pipeline ---
  dds   <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = t(matr), colData = meta, design = as.formula(form_str)))
  fit   <- suppressMessages(DESeq2::DESeq(dds, test = "Wald"))
  coefs <- suppressMessages(DESeq2::resultsNames(fit))


  # Normalized counts
  norm_counts  <- suppressMessages(counts(fit, normalized = TRUE))
  ave_expr_vec <- rowMeans(log2(norm_counts + 1))


  ## --- Contrast matrix again (DESeq internal) ---
  contrast_calculations <- list(
    baseline_1     = list(c(coefs[3]), c(coefs[1])),  # G1.A2 - G1.A1
    baseline_2     = list(c(coefs[4]), c(coefs[2])),  # G2.A2 - G2.A1
    relationship_1 = list(c(coefs[2]), c(coefs[1])),  # G2.A1 - G1.A1
    relationship_2 = list(c(coefs[4]), c(coefs[3])),  # G2.A2 - G1.A2
    interaction    = list(c(coefs[4], coefs[1]), c(coefs[3], coefs[2]))  
    # (G2.A2 - G1.A2) - (G2.A1 - G1.A1)
  )

  ## --- Extract results for each contrast ---
  stopifnot(
    length(names(contrast_calculations)) == length(names(contrast_calculations_pretty))
  )

  results_list <- lapply(names(contrast_calculations), function(cn) {
    contrast <- contrast_calculations[[cn]]

    res <- results(fit, contrast = contrast)
    tt  <- as.data.frame(res)


    data.frame(
      coef_id   = cn,
      coef_type = sub("_[0-9]+$", "", cn),
      contrast  = contrast_calculations_pretty[[cn]],
      g_1       = g_1,
      g_2       = g_2,
      a_1       = a_1,
      a_2       = a_2,
      feature   = rownames(tt),
      T_obs     = tt$log2FoldChange,
      SE        = tt$lfcSE,
      p_value   = tt$pvalue,
      p_adj     = tt$padj,
      ave_expr  = ave_expr_vec[rownames(tt)],
      row.names = NULL
    )
  })

  summary_stats <- do.call(rbind, results_list)

  return(summary_stats)
}
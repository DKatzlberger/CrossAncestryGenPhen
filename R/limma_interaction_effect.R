#' Limma-Based Interaction Test
#'
#' @param X Expression matrix for ancestry X. Rows = samples, columns = genes.
#' @param Y Expression matrix for ancestry Y. Rows = samples, columns = genes.
#' @param MX Metadata for X. Must include group and ancestry columns.
#' @param MY Metadata for Y. Must include group and ancestry columns.
#' @param g_col Name of the column indicating group (factor 1).
#' @param a_col Name of the column indicating ancestry (factor 2).
#' @param covariates Optional vector of covariate column names to adjust for.
#' @param verbose Logical, whether to print messages.
#'
#' @return summary stats
#' 
#' @export
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom limma voom lmFit eBayes topTable makeContrasts contrasts.fit
limma_interaction_effect <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col,
  covariates = NULL,
  use_voom = TRUE,
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
    .fun = "limma_interaction_effect"
  )


  ## --- Combine expression and metadata ---
  matr <- rbind(X, Y)
  meta <- rbind(MX, MY)

  if (!identical(rownames(matr), rownames(meta))) {
    stop("[limma_interaction_effect] Matrix and meta rownames must match exactly.")
  }


  ## --- Factor setup ---
  a_1 <- unique(MX[[a_col]])
  a_2 <- unique(MY[[a_col]])

  meta[[a_col]] <- factor(meta[[a_col]], levels = c(a_1, a_2))

  g_levels <- levels(meta[[g_col]])
  a_levels <- levels(meta[[a_col]])

  if (length(g_levels) != 2 || length(a_levels) != 2) {
    stop("[limma_interaction_effect] Function currently supports only 2x2 designs (two levels in g_col × two levels a_col).")
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

  ## --- Validation ---
  if (!identical(rownames(meta), rownames(design))) {
    matr <- matr[rownames(design), , drop = FALSE]
  }

  if (!identical(rownames(matr), rownames(design))) {
    stop("[limma_interaction_effect] Response matrix does not match design. Can happen if design matrix reorders sample values.")
  }


  ## --- Group vs covariates coef ---
  all_coefs   <- colnames(design)
  group_mask  <- grepl("^groups", all_coefs)
  group_coefs <- all_coefs[group_mask]
  covar_coefs <- all_coefs[!group_mask]


  ## --- Clean coef names ---
  group_coefs <- gsub("groups", "", group_coefs)               
  colnames(design)[group_mask] <- group_coefs

  # if (!is.null(covariates)) {
  #   for (cov in covariates) {
  #     if (is.factor(meta[[cov]]) || is.character(meta[[cov]])) {
  #       covar_coefs <- sub(paste0("^", cov), "", covar_coefs)
  #     }
  #   }
  #   colnames(design)[!group_mask] <- covar_coefs
  # }


  ## --- Define contrasts ---
  cols <- colnames(design)[group_mask]

  # Contrast calculations
  contrast_calculations <- list(
    baseline_1     = paste(cols[3], "-", cols[1]), # G1.A2 - G1.A1
    baseline_2     = paste(cols[4], "-", cols[2]), # G2.A2 - G2.A1
    relationship_1 = paste(cols[2], "-", cols[1]), # G2.A1 - G1.A1
    relationship_2 = paste(cols[4], "-", cols[3]), # G2.A2 - G1.A2
    interaction    = paste0("(", cols[4], " - ", cols[3], ") - (", cols[2], " - ", cols[1], ")")
  )

  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_calculations,
    levels    = design
  )


  ## --- Verbose message ---
  if (verbose) {
    message("\nLinear model summary:")
    message(sprintf("%-20s  %s", "Formula:", form_str))
    message(sprintf("%-20s  %s", "Groups:", paste(group_coefs, collapse = "  ")))
    message(sprintf("%-20s  %s", "Baseline:", paste(contrast_calculations[1:2], collapse = "  ")))
    message(sprintf("%-20s  %s", "Relationship:", paste(contrast_calculations[3:4], collapse = "  ")))
    message(sprintf("%-20s  %s", "Interaction:", paste(contrast_calculations[[5]], collapse = "  ")))
  }


  ## --- Model fit with/without voom ----
  if (use_voom) {
    dge <- edgeR::DGEList(counts = t(matr)) # samples x genes -> genes x samples
    dge <- edgeR::calcNormFactors(dge) 
    v   <- limma::voom(dge, design = design, plot = FALSE)
    fit <- limma::lmFit(v, design)
  } else {
    fit <- limma::lmFit(t(matr), design)                 
  }


  ## --- Apply contrasts + eBayes (always) ---
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)


  ## --- Extract results ---
  results_list <- lapply(
      seq_along(colnames(contrast_matrix)), function(i) {
      cn <- colnames(contrast_matrix)[i]
      tt <- limma::topTable(fit2, coef = cn, number = Inf, sort.by = "none")
      # Compute SE from logFC and t-stat
      tt$SE <- abs(tt$logFC / tt$t)

      data.frame(
        coef_id   = names(contrast_calculations)[i],
        coef_type = sub("_[0-9]+$", "", names(contrast_calculations)[i]),
        contrast  = cn,
        g_1       = g_1,
        g_2       = g_2,
        a_1       = a_1,
        a_2       = a_2,
        feature   = rownames(tt),
        T_obs     = tt$logFC,
        SE        = tt$SE,
        p_value   = tt$P.Value,
        p_adj     = tt$adj.P.Val,
        ave_expr  = tt$AveExpr,
        row.names = NULL
      )
    }
  )
  summary_stats <- do.call(rbind, results_list)


  ## --- Return ----
  return(summary_stats)
}

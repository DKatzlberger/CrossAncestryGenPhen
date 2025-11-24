#' Limma-Based Relationship Test on three datasets.
#'
#' @param R Expression matrix for refeernce ancestry. Rows = samples, columns = genes.
#' @param X Expression matrix for ancestry X. Rows = samples, columns = genes.
#' @param Y Expression matrix for ancestry Y. Rows = samples, columns = genes.
#' @param MR Metadata for R. Must include group and ancestry columns.
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
limma_correlation_effect <- function(
  R,
  X,
  Y,
  MR,
  MX,
  MY,
  g_col,
  a_col,
  covariates = NULL,
  use_voom = TRUE,
  verbose = TRUE
){

  ## --- Input data structure check ---
  # assert_input(
  #   R = R,
  #   X = X, 
  #   Y = Y,
  #   MR = MR,
  #   MX = MX, 
  #   MY = MY,
  #   g_col = g_col, 
  #   a_col = a_col,
  #   .fun = "limma_correlation_effect"
  # )


  ## --- Check data leakage ---
  if (length(intersect(rownames(R), rownames(X))) > 0) stop("Data leakage: R and X share rownames.")
  if (length(intersect(rownames(R), rownames(Y))) > 0) stop("Data leakage: R and Y share rownames.")
  if (length(intersect(rownames(X), rownames(Y))) > 0) stop("Data leakage: X and Y share rownames.")


  ## --- Ancestry validation ---
  a_1   <- unique(MX[[a_col]])
  a_2   <- unique(MY[[a_col]])
  a_1_1 <- unique(MR[[a_col]])
  if (a_1 != a_1_1) stop("[limma_correlation_effect] Ancestry level must be the same in Reference (R) as in Subset (X).")


  ## --- Combine expression and metadata ---
  expr_list <- list(R = R, X = X, Y = Y)
  meta_list <- list(R = MR, X = MX, Y = MY)

  ## --- Build the recipe ---
  form_str <- paste("~0 + groups")
  if (!is.null(covariates)) {
    form_str <- paste(form_str, "+", paste(covariates, collapse = " + "))
  }


  ## --- Relationship loop ---
  results_list <- lapply(names(expr_list), function(a_name) {

    matr <- expr_list[[a_name]]
    meta <- meta_list[[a_name]]

    if (!identical(rownames(matr), rownames(meta))) {
      stop(sprintf("[limma_correlation_effect] Matrix and meta rownames must match exactly."))
    }


    ## --- Ensure 2-level group ----
    if (!is.factor(meta[[g_col]])) stop("[limma_correlation_effect] g_col is not a factor.")
    g_levels <- levels(meta[[g_col]])
    a_levels <- unique(meta[[a_col]])

    if (length(g_levels) != 2) stop(sprintf("[limma_correlation_effect] Function currently supports only 2 groups (two levels in g_col)."))
    if (length(a_levels) != 1) stop(sprintf("[limma_correlation_effect] Function currently supports only 1 ancestry (one level in a_col)."))

    g_1 <- g_levels[1]
    g_2 <- g_levels[2]
    a_1 <- a_levels[1]

    meta[["groups"]] <- factor(
      paste(
        meta[[g_col]], 
        meta[[a_col]], 
        sep = "."
      ), 
      levels = c(
        paste(g_1, a_1, sep = "."),  
        paste(g_2, a_1, sep = ".")
      )
    )


    ## --- Means model ---
    design <- model.matrix(as.formula(form_str), data = meta)
    colnames(design) <- make.names(colnames(design))


    ## --- Validation ---
    if (!identical(rownames(meta), rownames(design))) {
      matr <- matr[rownames(design), , drop = FALSE]
    }

    if (!identical(rownames(matr), rownames(design))) {
      stop("[limma_correlation_effect] Response matrix does not match design. Can happen if design matrix reorders sample values.")
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

    contrast_calculations <- list(relationship = paste(cols[2], "-", cols[1]))
    contrast_matrix <- limma::makeContrasts(
      contrasts = contrast_calculations,
      levels = design
    )


    ## --- Model fit with/without voom ----
    if (use_voom) {
      dge <- edgeR::DGEList(counts = t(matr))
      dge <- edgeR::calcNormFactors(dge)
      v   <- limma::voom(dge, design, plot = FALSE)
      fit <- limma::lmFit(v, design)
    } else {
      fit <- limma::lmFit(t(matr), design)
    }


    ## --- Apply contrasts + eBayes (always) ---
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    fit2 <- limma::eBayes(fit2)


    ## --- Extract results ---
    tt <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")
    tt$SE <- abs(tt$logFC / tt$t)


    ## --- Return ---
    list(
      stats = data.frame(
        coef_id   = paste0("relationship_", a_name),
        coef_type = "relationship",
        contrast  = paste(cols[2], "-", cols[1]),
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
      ),
      coefs = group_coefs
    )
  })

  ## --- Merge all ancestry results ---
  summary_stats <- do.call(rbind, lapply(results_list, `[[`, "stats"))
  R_coefs <- results_list[[1]]$coefs
  X_coefs <- results_list[[2]]$coefs
  Y_coefs <- results_list[[3]]$coefs


  ## --- Verbose message ---
  if (verbose) {
    message("\nLinear model summary:")
    message(sprintf("%-20s  %s", "Formula:", form_str))
    message(sprintf("%-20s  %-30s Contrast: %s", paste0("Reference R (", a_1, "):"), paste(R_coefs[1], R_coefs[2]), paste(R_coefs[2], "-", R_coefs[1])))
    message(sprintf("%-20s  %-30s Contrast: %s", paste0("Subset    X (", a_1, "):"), paste(X_coefs[1], X_coefs[2]), paste(X_coefs[2], "-", X_coefs[1])))
    message(sprintf("%-20s  %-30s Contrast: %s", paste0("Inference Y (", a_2, "):"), paste(Y_coefs[1], Y_coefs[2]), paste(Y_coefs[2], "-", Y_coefs[1])))
  }

  ## --- Correlation ---
  compute_correlations <- function(summary_stats) {

    # split into R, X, Y
    R <- summary_stats[grepl("R$", summary_stats$coef_id), c("feature", "T_obs")]
    X <- summary_stats[grepl("X$", summary_stats$coef_id), c("feature", "T_obs")]
    Y <- summary_stats[grepl("Y$", summary_stats$coef_id), c("feature", "T_obs")]

    # merge on feature
    RX <- merge(R, X, by = "feature", suffixes = c("_R", "_X"))
    RY <- merge(R, Y, by = "feature", suffixes = c("_R", "_Y"))

    # define correlation methods
    cor_methods <- list(
      pearson  = function(a, b) cor(a, b, use = "pairwise", method = "pearson"),
      spearman = function(a, b) cor(a, b, use = "pairwise", method = "spearman")
    )

    # compute correlations
    out <- lapply(names(cor_methods), function(m) {
      fun <- cor_methods[[m]]

      cor_RX <- fun(RX$T_obs_R, RX$T_obs_X)
      cor_RY <- fun(RY$T_obs_R, RY$T_obs_Y)

      list(
        method = m,
        cor_RX = cor_RX,
        cor_RY = cor_RY
      )
    })

    return(out)
  }
  cors <- compute_correlations(summary_stats)

  if (verbose){
    message("\nCorrelation performance:")
    message(sprintf("%-20s  Pearson %.4f   Spearman %.4f ", paste0("Subset    X (", a_1, "):"), cors[[1]]$cor_RX, cors[[2]]$cor_RX))
    message(sprintf("%-20s  Pearson %.4f   Spearman %.4f ", paste0("Inference Y (", a_2, "):"), cors[[1]]$cor_RY, cors[[2]]$cor_RY))
  }


  ## --- Return ---
  return(summary_stats)
}

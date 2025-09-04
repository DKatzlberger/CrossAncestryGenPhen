#' Internal: minimal checks for X/Y/R (+ optional MX/MY/MR, g_col, a_col)
#' Skips any argument left as NULL.
#' @noRd
assert_input <- function(
  X  = NULL, 
  Y  = NULL, 
  R  = NULL,
  MX = NULL, 
  MY = NULL,
  MR = NULL,
  g_col = NULL, 
  a_col = NULL
) {
  ## Matrices + rownames
  if (!is.null(X)) {
    if (!is.matrix(X)) stop("X must be a matrix.")
    if (is.null(rownames(X))) stop("X must have rownames corresponding to sample IDs.")
  }
  if (!is.null(Y)) {
    if (!is.matrix(Y)) stop("Y must be a matrix.")
    if (is.null(rownames(Y))) stop("Y must have rownames corresponding to sample IDs.")
  }
  if (!is.null(R)) {
    if (!is.matrix(R)) stop("R must be a matrix.")
    if (is.null(rownames(R))) stop("R must have rownames corresponding to sample IDs.")
  }

  ## Metadata frames + rownames
  if (!is.null(MX)) {
    if (!is.data.frame(MX)) stop("MX must be a data.frame.")
    if (is.null(rownames(MX))) stop("MX must have rownames corresponding to sample IDs.")
  }
  if (!is.null(MY)) {
    if (!is.data.frame(MY)) stop("MY must be a data.frame.")
    if (is.null(rownames(MY))) stop("MY must have rownames corresponding to sample IDs.")
  }
  if (!is.null(MR)) {
    if (!is.data.frame(MR)) stop("MR must be a data.frame.")
    if (is.null(rownames(MR))) stop("MR must have rownames corresponding to sample IDs.")
  }

  ## Row count alignment
  if (!is.null(X) && !is.null(MX) && nrow(X) != nrow(MX)) {
    stop("X/MX row counts must match.")
  }
  if (!is.null(Y) && !is.null(MY) && nrow(Y) != nrow(MY)) {
    stop("Y/MY row counts must match.")
  }
  if (!is.null(R) && !is.null(MR) && nrow(R) != nrow(MR)) {
    stop("R/MR row counts must match.")
  }

  ## Required columns
  md_list <- list(MX = MX, MY = MY, MR = MR)
  if (!is.null(g_col) || !is.null(a_col)) {
    for (nm in names(md_list)) {
      M <- md_list[[nm]]
      if (is.null(M)) next
      if (!is.null(g_col) && !(g_col %in% names(M))) {
        stop(sprintf("g_col must exist in %s.", nm))
      }
      if (!is.null(a_col) && !(a_col %in% names(M))) {
        stop(sprintf("a_col must exist in %s.", nm))
      }
    }
  }

  ## Ancestry uniqueness: per block only
  ancestry_labels <- NULL
  if (!is.null(a_col)) {
    labels <- list()
    for (nm in names(md_list)) {
      M <- md_list[[nm]]
      if (is.null(M)) next
      u <- unique(na.omit(M[[a_col]]))
      if (length(u) != 1L) {
        stop(sprintf("In %s, %s must have exactly one non-NA value.", nm, a_col))
      }
      labels[[nm]] <- as.character(u)
    }
    ancestry_labels <- labels
  }

  ## Group level consistency across all provided metadata
  if (!is.null(g_col)) {
    g_levels_list <- lapply(md_list[!vapply(md_list, is.null, logical(1))],
                            function(M) sort(unique(as.character(M[[g_col]]))))
    if (length(g_levels_list) > 1) {
      ref <- g_levels_list[[1]]
      for (nm in names(g_levels_list)[-1]) {
        if (!identical(g_levels_list[[nm]], ref)) {
          stop(sprintf("g_col levels differ between metadata frames. %s has %s, expected %s",
                       nm,
                       paste(g_levels_list[[nm]], collapse = ","),
                       paste(ref, collapse = ",")))
        }
      }
    }
  }

  invisible(ancestry_labels %||% TRUE)
}

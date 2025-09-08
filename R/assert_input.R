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
  a_col = NULL,
  .fun = NULL
) {
  .bold <- function(x) paste0("**", as.character(x), "**")
  .caller <- function() {
    if (!is.null(.fun)) return(as.character(.fun))
    cs <- sys.calls()
    if (length(cs) >= 2L) {
      as.character(cs[[length(cs) - 1L]][[1L]])
    } else {
      "assert_input"
    }
  }
  .stopf <- function(msg, ...) {
    args <- lapply(list(...), .bold)
    msg  <- do.call(sprintf, c(list(msg), args))
    stop(sprintf("[%s] %s", .caller(), msg), call. = FALSE)
  }

  ## Matrices + rownames
  if (!is.null(X)) {
    if (!is.matrix(X)) .stopf("X must be a matrix.")
    if (is.null(rownames(X))) .stopf("X must have rownames corresponding to sample IDs.")
  }
  if (!is.null(Y)) {
    if (!is.matrix(Y)) .stopf("Y must be a matrix.")
    if (is.null(rownames(Y))) .stopf("Y must have rownames corresponding to sample IDs.")
  }
  if (!is.null(R)) {
    if (!is.matrix(R)) .stopf("R must be a matrix.")
    if (is.null(rownames(R))) .stopf("R must have rownames corresponding to sample IDs.")
  }

  ## Metadata frames + rownames
  if (!is.null(MX)) {
    if (!is.data.frame(MX)) .stopf("MX must be a data.frame.")
    if (is.null(rownames(MX))) .stopf("MX must have rownames corresponding to sample IDs.")
  }
  if (!is.null(MY)) {
    if (!is.data.frame(MY)) .stopf("MY must be a data.frame.")
    if (is.null(rownames(MY))) .stopf("MY must have rownames corresponding to sample IDs.")
  }
  if (!is.null(MR)) {
    if (!is.data.frame(MR)) .stopf("MR must be a data.frame.")
    if (is.null(rownames(MR))) .stopf("MR must have rownames corresponding to sample IDs.")
  }

  ## Row count alignment
  if (!is.null(X) && !is.null(MX) && nrow(X) != nrow(MX)) {
    .stopf("X/MX row counts must match (got %s vs %s).", nrow(X), nrow(MX))
  }
  if (!is.null(Y) && !is.null(MY) && nrow(Y) != nrow(MY)) {
    .stopf("Y/MY row counts must match (got %s vs %s).", nrow(Y), nrow(MY))
  }
  if (!is.null(R) && !is.null(MR) && nrow(R) != nrow(MR)) {
    .stopf("R/MR row counts must match (got %s vs %s).", nrow(R), nrow(MR))
  }

  ## Required columns
  md_list <- list(MX = MX, MY = MY, MR = MR)
  if (!is.null(g_col) || !is.null(a_col)) {
    for (nm in names(md_list)) {
      M <- md_list[[nm]]
      if (is.null(M)) next
      if (!is.null(g_col) && !(g_col %in% names(M))) {
        .stopf("g_col must exist in %s (missing %s).", nm, g_col)
      }
      if (!is.null(a_col) && !(a_col %in% names(M))) {
        .stopf("a_col must exist in %s (missing %s).", nm, a_col)
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
        .stopf("In %s, %s must have exactly one non-NA value (got %s).",
               nm, a_col, paste(u, collapse = ", "))
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
          .stopf("g_col levels differ between metadata frames. %s has %s, expected %s",
                 nm,
                 paste(g_levels_list[[nm]], collapse = ","),
                 paste(ref, collapse = ","))
        }
      }
    }
  }

  invisible(ancestry_labels %||% TRUE)
}

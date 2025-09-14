#' Build a parameter grid
#'
#' Creates a grid of parameter combinations with `expand.grid()`.  
#' Optionally enforces a rule for `n_degs` and `log2fc` (set both to 0 if either is 0),  
#' and prints the grid when `verbose = TRUE`.
#'
#' @param ... Named vectors of parameters.
#' @param zero_rules Logical, default `TRUE`. Apply `n_degs`/`log2fc` zeroing rule.
#' @param verbose Logical, default `TRUE`. Print the grid and number of rows.
#'
#' @return A data.frame of parameter combinations.
#'
#' @examples
#' sim_make_grid(n_samples = c(10), n_degs = c(0,1000), log2fc = c(0,1))
#' sim_make_grid(total_samples = 1000, between_ratio = 2, within_ratio = 2, zero_rules = FALSE)
#'
#' @export
sim_make_grid <- function(
  ..., 
  zero_rules = TRUE,
  verbose = TRUE
) {

  ## --- Collect arguments in a list ---
  args <- list(...)
  

  ## --- Expand the grid ---
  grid <- expand.grid(args)
  

  ## --- Remove the zero rule ---
  if (zero_rules) {
    if (all(c("n_degs", "log2fc") %in% names(grid))) {
      zero_idx <- grid$n_degs == 0 | grid$log2fc == 0
      grid$n_degs[zero_idx] <- 0
      grid$log2fc[zero_idx] <- 0
      grid <- unique(grid)
    }
  }
  
  ## --- Verbose message ---
  if (verbose){
    message("\nGrid summary:")
    message("Unique combinations: ", nrow(grid))
    msg <- capture.output(print(grid, row.names = FALSE))
    message(paste(msg, collapse = "\n"))
  }
  
  ## --- Return ---
  return(grid)
}

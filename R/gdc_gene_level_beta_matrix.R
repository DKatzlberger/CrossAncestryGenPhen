#' Create GDC gene-level beta matrix
#'
#' Filters a methylation matrix to probes present in a
#' probe-to-gene mapping table and removes probes with
#' any missing values. Returns a cleaned beta matrix.
#'
#' @param meth Input beta matrix (samples × probes)
#' @param probe2gene Probe-to-gene mapping table
#' @param verbose Print progress messages
#'
#' @return Filtered beta matrix
#'
#' @importFrom stats na.omit
#' @export
#'
gdc_gene_level_beta_matrix <- function(
  meth,
  probe2gene,
  verbose = TRUE       
){

  ## --- Extract probes ---
  probes_meth <- colnames(meth)
  probes_map  <- probe2gene$probe


  ## --- Filter matrix ---
  keep <- intersect(probes_meth, probes_map)
  meth <- meth[, keep, drop = FALSE]
  if (verbose) message("Data: ", length(probes_meth), " | Map: ", length(probes_map), " | Overlap: ", length(keep))


  ## --- Remove probes with ANY NA values ---
  isNA <- colSums(is.na(meth)) != 0
  meth <- meth[, !isNA, drop = FALSE]
  if (verbose) if (verbose) message("Removed ", sum(isNA), " probes with NA values")


  ## --- Verbose message ---
  if (verbose){
    message("\nGDC gene-level beta matrix summary:")
    message(sprintf("Dimensions: %d (samples) x %d (features)", nrow(meth), ncol(meth)))
  }


  ## --- Return ---
  return(meth)
}
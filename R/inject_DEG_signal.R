#' Inject Simulated DEGs and Interaction Effects into NB Counts
#'
#' Takes output from `simulate_NB_counts()` and injects differential expression
#' by adjusting counts for selected genes and groups.
#'
#' @param sim_data Output from `simulate_NB_counts()`.
#' @param n_degs Total number of DEGs to inject.
#' @param log2fc Log2 fold change for DEGs.
#' @param prop_interaction Proportion of DEGs that should have interaction effects.
#' @param seed Optional random seed.
#'
#' @return A list with updated counts, sample info, and DEG metadata.
#' @export
inject_DEG_signal <- function(
  sim_data,
  n_degs = 1000,
  log2fc = 1.5,
  prop_interaction = 0.3,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  counts <- sim_data$counts  # genes × samples
  info <- sim_data$sample_info
  params <- sim_data$params

  all_genes <- rownames(counts)
  n_interaction <- round(n_degs * prop_interaction)
  n_main_each <- round((n_degs - n_interaction) / 2)

  deg_genes <- sample(all_genes, n_degs)
  main_a1_genes <- deg_genes[1:n_main_each]
  main_a2_genes <- deg_genes[(n_main_each + 1):(2 * n_main_each)]
  interaction_genes <- deg_genes[(2 * n_main_each + 1):n_degs]

  # Convenience
  fc <- 2^log2fc

  # Helper to get column indices (samples) by group
  group_cols <- function(g) {
    which(info$group == g)
  }

  # Inject main effects
  inject_fc <- function(genes, g2, flip = FALSE) {
    for (gene in genes) {
      gene_idx <- which(rownames(counts) == gene)
      sample_idx <- group_cols(g2)

      mu <- params$mu[params$feature == gene]
      disp <- params$dispersion[params$feature == gene]
      mu_fc <- mu * ifelse(flip, 1 / fc, fc)

      counts[gene_idx, sample_idx] <- rnbinom(length(sample_idx), mu = mu_fc, size = disp)
    }
  }

  inject_fc(main_a1_genes, g2 = "a1_g2", flip = FALSE)
  inject_fc(main_a2_genes, g2 = "a2_g2", flip = FALSE)

  # Inject interaction: opposite effects for a1_g2 and a2_g2
  for (gene in interaction_genes) {
    gene_idx <- which(rownames(counts) == gene)
    mu <- params$mu[params$feature == gene]
    disp <- params$dispersion[params$feature == gene]

    idx_a1 <- group_cols("a1_g2")
    idx_a2 <- group_cols("a2_g2")

    counts[gene_idx, idx_a1] <- rnbinom(length(idx_a1), mu = mu * fc, size = disp)
    counts[gene_idx, idx_a2] <- rnbinom(length(idx_a2), mu = mu / fc, size = disp)
  }

  return(list(
    counts = counts,
    sample_info = info,
    degs = list(
      main_a1 = main_a1_genes,
      main_a2 = main_a2_genes,
      interaction = interaction_genes,
      log2fc = log2fc
    )
  ))
}

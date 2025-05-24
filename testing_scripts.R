library(CrossAncestryGenPhen)
library(ggplot2)

set.seed(1)
n_per_group <- 25
p <- 2000

# Simulate expression data
X <- matrix(rnorm(n_per_group * 2 * p), nrow = n_per_group * 2, ncol = p)
Y <- matrix(rnorm(n_per_group * 2 * p), nrow = n_per_group * 2, ncol = p)
colnames(X) <- colnames(Y) <- paste0("Gene_", seq_len(p))

# Simulate metadata
MX <- data.frame(
condition = factor(rep(c("A", "B"), each = n_per_group)),
ancestry = "EUR"
)
MY <- data.frame(
condition = factor(rep(c("A", "B"), each = n_per_group)),
ancestry = "AFR"
)

# Run the function from your installed package
result <- CrossAncestryGenPhen::perm_diff_interaction(
  X = X,
  Y = Y,
  MX = MX,
  MY = MY,
  g_col = "condition",
  a_col = "ancestry",
  B = 100,
  seed = 42
)

p <- plot_perm_distribution(
  result = result,
  features = c("Gene_1", "Gene_2", "Gene_3", "Gene_4", "Gene_5")
  )

ggsave("perm_plot.pdf", p, width = 8, height = 4)
# Check if the plot was saved correctly



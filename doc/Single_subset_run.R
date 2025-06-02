## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data setup, fig.alt = "Bar plot showing condition imbalance across ancestries"----
library(CrossAncestryGenPhen)
library(ggplot2)

# Seed for reproducibility
set.seed(42)

# Simulate example data
p <- 100  # Number of genes
n_EUR <- 600 
n_AFR <- 40  

# Expression matrices for EUR and AFR ancestries
X <- matrix(rnorm(n_EUR * p), nrow = n_EUR, ncol = p)
Y <- matrix(rnorm(n_AFR * p), nrow = n_AFR, ncol = p)
colnames(X) <- colnames(Y) <- paste0("Gene_", seq_len(p))

# Metadata for EUR and AFR ancestries
# EUR: overrepresented compared to AFR
MX <- data.frame(
  condition = factor(c(rep("H", 400), rep("D", 200))),
  ancestry = "EUR"
)

# AFR: underrepresented compared to EUR
MY <- data.frame(
  condition = factor(c(rep("H", 10), rep("D", 30))),
  ancestry = "AFR"
)

# Visualize sample size imbalance
meta <- rbind(MX, MY)

# Plot
ggplot(meta, aes(x = ancestry, fill = condition)) +
  geom_bar(position = "dodge", color = "black") +
  labs(
    title = "Condition Imbalance Across Ancestries",
    x = "Ancestry",
    y = "Sample Count",
    fill = "Condition"
  ) 

## ----stratification, fig.alt = "Bar plot showing balanced subset"-------------
# Define column to stratify
stratify_cols <- c("condition")

# Split the data into stratified sets
split <- split_stratified_ancestry_sets(
  X = X,
  Y = Y,
  MX = MX,
  MY = MY,
  stratify_cols = stratify_cols,
  seed = 42
)

# Visulaize stratified sets
plot_stratified_sets(
  x = split,
  stratify_cols = stratify_cols
)

## ----split_output-------------------------------------------------------------
# Output
str(split)

## ----interaction_analysis, fig.alt = "Bar plot showing interaction results"----
# Define the ancestry column and the condition column
g_col <- "condition"
a_col <- "ancestry"

# Define the number of permutations
B <- 1000

result <- perm_diff_interaction(
  X = split$test$X,
  Y = split$inference$X,
  MX = split$test$M,
  MY = split$inference$M,
  g_col = g_col,
  a_col = a_col,
  permute = TRUE,   # Perform permutation
  B = B,            # Number of iterations (shared for both bootstrapping and permutation)
  seed = 42,
  check_convergence = TRUE
)

## ----summary_stats------------------------------------------------------------
head(result$summary_stats, 10)

## ----B_used-------------------------------------------------------------------
result$B_used_boot
result$B_used_perm

## ----T_boot_distribution, fig.alt = "Distribution of bootstrapped test statistics"----
plot_T_distribution(
  x = result, 
  statistic = "T_boot",
  title = "Bootstrapped distribution with 95% confidence intervals"
)

## ----T_perm_distribution, fig.alt = "Distribution of permuted test statistics"----
plot_T_distribution(
  x = result, 
  statistic = "T_perm",
  title = "Permuted distribution with two-sided p-values"
)


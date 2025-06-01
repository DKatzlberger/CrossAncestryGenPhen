## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data setup---------------------------------------------------------------
library(CrossAncestryGenPhen)

# Seed for reproducibility
set.seed(42)

# Simulate example data
p <- 1000 # Number of genes
n_EUR <- 120 # Number of individuals in EUR ancestry
n_AFR <- 40  # Number of individuals in AFR ancestry

# Expression matrices for EUR and AFR ancestries
X <- matrix(rnorm(n_EUR * p), nrow = n_EUR, ncol = p)
Y <- matrix(rnorm(n_AFR * p), nrow = n_AFR, ncol = p)
colnames(X) <- colnames(Y) <- paste0("Gene_", seq_len(p))

# Metadata for EUR and AFR ancestries
# EUR: overrepresented compared to AFR
MX <- data.frame(
  condition = factor(c(rep("H", 70), rep("D", 30))),
  ancestry = "EUR"
)

# AFR: underrepresented compared to EUR
MY <- data.frame(
  condition = factor(c(rep("H", 10), rep("D", 30))),
  ancestry = "AFR"
)

# Visualize sample size imbalance
library(ggplot2)

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

## ----stratification-----------------------------------------------------------



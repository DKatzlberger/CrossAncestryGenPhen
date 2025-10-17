#' Plot sample expression density
#'
#' Plot density curves of expression values for selected samples.
#'
#' @param X Numeric matrix or data frame for the first group (samples in rows, genes in columns).
#' @param Y Numeric matrix or data frame for the second group (same orientation as X).
#' @param samples Vector of sample indices or row names to plot. If NULL, the first 9 samples are used.
#' @param cpm Logical; if TRUE, transform counts to log2 CPM before plotting.
#' @param mval Logical; if TRUE, apply m-value transformation before t-SNE.
#' @param title Plot title.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom edgeR cpm
#' @importFrom ggplot2 ggplot aes geom_density labs
plot_sample_density <- function(
  X,
  Y,
  samples = NULL,
  cpm = FALSE,
  mval = FALSE,
  title = NULL,
  x_label = NULL,
  y_label = NULL 
){

  ## --- Input data structure check ---
  assert_input(
    X = X,
    Y = Y,
    .fun = "plot_sample_density"
  )


  ## --- Combine matrices ---
  X_comb <- rbind(X, Y)


  ## --- CPM transform if requested ---
  if (cpm) {
    X_cpm <- edgeR::cpm(t(X_comb), log = TRUE)
    X_comb <- t(X_cpm)  
  }

  if (mval){
    X_mval <- beta_to_mval(t(X_comb))
    X_comb <- t(X_mval)
  }


  ## --- Select samples ---
  if (is.null(samples)) {
    sample_X <- head(rownames(X), 5)
    sample_Y <- head(rownames(Y), 5)
    samples  <- c(sample_X, sample_Y)
  } else {
    if (is.character(samples)) {
      samples <- match(samples, rownames(X_comb))
    }
  }

  X_sel <- X_comb[samples, , drop = FALSE]
  sel_names <- rownames(X_sel)


  ## --- Reshape with base R (stack) ---
  df_long <- stack(as.data.frame(t(X_sel)))  
  df_long$Sample <- rep(sel_names, each = ncol(X_sel))
  colnames(df_long) <- c("Value", "Gene", "Sample")


  ## --- Plot ---
  p <- ggplot(
    data = df_long, 
    mapping = aes(
      x = Value, 
      color = Sample
      )
    ) +
    geom_density() +
    labs(
      title = title,
       x = ifelse(!is.null(x_label), x_label,if (cpm) {"Log2 CPM"} else if (mval) {"M-values"} else {"Raw values"}),
      y = ifelse(is.null(y_label), "Density", y_label),
      color = "Sample"
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  ## --- Return ---
  return(p)
}
#' Plot concordance of p-values between two methods
#'
#' Merge two result tables by feature and compare their p-values (or other
#' numeric columns) on a \eqn{-\log_{10}} scale using a scatter or hexbin plot.
#' Adds a 1:1 reference line and a fitted regression line. Optionally label
#' selected features.
#'
#' @param x_data Data frame with a `feature` column and a numeric column given in `x_var`.
#' @param y_data Data frame with a `feature` column and a numeric column given in `y_var`.
#' @param x_var Character, name of the p-value column in `x_data`.
#' @param y_var Character, name of the p-value column in `y_data`.
#' @param features Optional character vector of feature IDs to label.
#' @param hex Logical; if `TRUE` (default) draw a hexbin plot, otherwise scatter points.
#' @param title Optional plot title.
#' @param x_label,y_label Axis labels (defaults use `x_var`/`y_var`).
#' @param point_size Numeric; point size if `hex = FALSE`. Default is 1.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_pvalue_concordance(
#'   x_data = edgeR_res[edgeR_res$coef_id == "interaction", ],
#'   y_data = limma_res[limma_res$coef_id == "interaction", ],
#'   x_var = "p_value",
#'   y_var = "p_value",
#'   title = "Interaction concordance",
#'   features = c("TP53", "BRCA1"),
#'   hex = FALSE
#' )
#' }
plot_pvalue_concordance <- function(
  x_data,
  y_data,
  x_var,
  y_var,
  features = NULL,
  hex = TRUE,
  title = NULL, 
  x_label = NULL,
  y_label = NULL,
  point_size = 1
){
  x_data <- x_data[, c("feature", x_var)]
  y_data <- y_data[, c("feature", y_var)]

  ## --- Make wide ---
  wide_data <- merge(
    x = x_data, 
    y = y_data,
    by = c("feature"),
    suffixes = c("_x", "_y")
  )

  ## --- Annotate features if provided and 'feature' exists ---
  if (!is.null(features) && "feature" %in% colnames(wide_data)) {
    wide_data$annotate <- ifelse(wide_data$feature %in% features, wide_data$feature, NA)
  } else {
    wide_data$annotate <- NA
  }

  ## --- Make symbol ---
  x_col <- sym(paste0(x_var, "_x"))
  y_col <- sym(paste0(y_var, "_y"))

  ## --- Base plot ---
  p <- ggplot(
    data = wide_data,
    mapping = aes(
      x = !!x_col,
      y = !!y_col
    )
  ) +
  geom_abline(
    slope = 1,
    linewidth = 0.3,
    color = "grey50"
  )

  ## --- Decide on geom ---
  if (hex) {
    p <- p + geom_hex(
      linewidth = point_size
    )
  } else {
    p <- p + geom_point(
      size = point_size
    )
  }

  ## --- Add correlation line ---
  p <- p + geom_smooth(
    method = "lm",
    se = FALSE, 
    color = "blue",
    linewidth = 0.3
  )

  ## --- Axis scale ---
  neglog10_trans <- scales::trans_new(
    name      = "neglog10",
    transform = function(x) -log10(x),
    inverse   = function(x) 10^(-x),
    domain    = c(1e-300, 1),
    breaks    = scales::trans_breaks(function(x) -log10(x), function(x) 10^(-x)),
    format    = scales::scientific_format(digits = 1)
  )

  p <- p + scale_x_continuous(
    trans = neglog10_trans
  ) +
  scale_y_continuous(
    trans = neglog10_trans
  )


  ## --- Feature annotation ---
  if (!is.null(features) && "feature" %in% colnames(df)) {
    p <- p + ggrepel::geom_text_repel(
      data = subset(df, !is.na(annotate)),
      aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
        label = annotate
      ),
      size = 2,
      min.segment.length = 0,
      segment.size = 0.3,
      max.overlaps = Inf,
      force = 2,
      box.padding = 0.5,
      bg.color = "white",
      bg.r = 0.2,
      inherit.aes = FALSE
    )
  }

  ## --- Final styling ---
  p <- p +
  labs(
    title = title,
    x = ifelse(is.null(x_label), x_var, x_label),
    y = ifelse(is.null(y_label), y_var , y_label)
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_small_legend() 


  return(p)
}
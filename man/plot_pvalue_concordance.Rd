% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pvalue_concordance.R
\name{plot_pvalue_concordance}
\alias{plot_pvalue_concordance}
\title{Plot concordance of -log10 p-values between two methods}
\usage{
plot_pvalue_concordance(
  data,
  x_var,
  y_var,
  x_sig_source = NULL,
  y_sig_source = NULL,
  facet_col = NULL,
  facet_levels = NULL,
  x_label = paste("-log10(", x_var, ")"),
  y_label = paste("-log10(", y_var, ")"),
  title = NULL,
  sig_thr = 0.05,
  epsilon = 1e-16,
  log_cap = 5,
  point_size = 0.5,
  hex = FALSE
)
}
\arguments{
\item{data}{Data frame containing p-values and optional grouping or significance columns.}

\item{x_var}{Character. Name of the column for x-axis p-values.}

\item{y_var}{Character. Name of the column for y-axis p-values.}

\item{x_sig_source}{Optional character. Column for significance of x p-values.}

\item{y_sig_source}{Optional character. Column for significance of y p-values.}

\item{facet_col}{Optional character. Column name to facet the plot by.}

\item{facet_levels}{Optional character vector. Subset of levels to include from \code{facet_col}.}

\item{x_label}{Character. X-axis label. Defaults to \code{-log10(x_var)}.}

\item{y_label}{Character. Y-axis label. Defaults to \code{-log10(y_var)}.}

\item{title}{Character. Optional plot title.}

\item{sig_thr}{Numeric. P-value threshold for calling significance. Default is 0.05.}

\item{epsilon}{Numeric. Small constant added before log10 transform. Default is 1e-16.}

\item{log_cap}{Numeric. Maximum value shown for -log10(p). Default is 5.}

\item{point_size}{Numeric. Size of points in the scatter plot. Default is 0.5.}

\item{hex}{Logical. If TRUE, use hexagonal binning instead of individual points. Default is FALSE.}
}
\value{
A ggplot2 object showing the concordance of p-values between two methods.
}
\description{
Creates a scatter or hexbin plot comparing -log10 transformed p-values from two
statistical methods. Optionally colors points by significance and facets
by a grouping variable. Correlation coefficients and regression lines are shown
globally or per facet.
}
\details{
If \code{hex = TRUE}, a hexbin plot is produced without individual points or significance coloring.
Otherwise, a scatter plot is shown with optional significance-based coloring and shape encoding for capped values.
Correlation (Pearson's r) and regression lines are shown globally or per facet if faceting is enabled.
}

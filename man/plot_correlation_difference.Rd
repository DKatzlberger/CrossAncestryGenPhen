% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_correlation_difference.R
\name{plot_correlation_difference}
\alias{plot_correlation_difference}
\title{Plot correlation differences with optional facets}
\usage{
plot_correlation_difference(
  x,
  cor_cols = c("XR", "YR"),
  t_col,
  p_col,
  facet_col = NULL,
  facet_levels = NULL,
  x_label = NULL,
  y_label = NULL,
  title = NULL,
  point_size = 0.5
)
}
\arguments{
\item{x}{A data.frame or data.table with required columns.}

\item{cor_cols}{Two columns with correlation values to compare.}

\item{t_col}{Column with T statistic values for annotation.}

\item{p_col}{Column with p-values for annotation.}

\item{facet_col}{Optional column to use for faceting.}

\item{facet_levels}{Optional facet values to include in plot.}

\item{x_label}{Label for the x-axis (default is NULL).}

\item{y_label}{Label for the y-axis (default is NULL).}

\item{title}{Optional plot title.}
}
\value{
A ggplot2 object showing the correlation comparison.
}
\description{
Create a bar plot comparing correlation values across conditions.
}

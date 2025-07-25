% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot_correlation_diff.R
\name{boot_correlation_diff}
\alias{boot_correlation_diff}
\title{Bootstrap-Based Correlation Difference Test (Unstratified)}
\usage{
boot_correlation_diff(
  X,
  Y,
  R,
  MX,
  MY,
  MR,
  g_col,
  method = c("pearson", "spearman"),
  B = 1000,
  seed = NULL,
  alpha = 0.05
)
}
\arguments{
\item{X}{Expression matrix for ancestry X (samples x genes)}

\item{Y}{Expression matrix for ancestry Y}

\item{R}{Expression matrix for reference ancestry}

\item{MX}{Metadata for X}

\item{MY}{Metadata for Y}

\item{MR}{Metadata for R}

\item{g_col}{Column in metadata specifying condition/group (factor with 2 levels)}

\item{method}{Correlation method: "pearson" (default) or "spearman"}

\item{B}{Number of bootstrap samples}

\item{seed}{Optional seed for reproducibility}

\item{alpha}{Significance level (e.g., 0.05 for 95\% CI)}
}
\description{
Estimates uncertainty in the correlation difference between ancestries X and Y
(relative to reference R), using non-stratified bootstrap resampling within ancestries.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/track_sample_ids.R
\name{track_sample_ids}
\alias{track_sample_ids}
\title{Track Sample Roles and IDs from a Split}
\usage{
track_sample_ids(split, iteration)
}
\arguments{
\item{split}{A split object returned from \code{split_stratified_ancestry_sets()}.}

\item{iteration}{An integer iteration number to annotate the split with.}
}
\value{
A data.frame with columns: \code{ids}, \code{role}, and \code{iteration}.
}
\description{
Extracts the sample IDs from a split object and labels their roles (test, train, inference),
associating each sample with its assigned role in a given iteration.
}

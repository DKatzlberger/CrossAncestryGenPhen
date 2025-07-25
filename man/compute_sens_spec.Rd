% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_spec_sens.R
\name{compute_sens_spec}
\alias{compute_sens_spec}
\title{Compute sensitivity and specificity for feature selection}
\usage{
compute_sens_spec(true, selected, all)
}
\arguments{
\item{true}{Character vector of true significant features (ground truth).}

\item{selected}{Character vector of features selected by your method.}

\item{all}{Character vector of all features tested or eligible to select.}
}
\value{
A data frame with TP, FP, FN, TN, sensitivity, and specificity.
}
\description{
Compares selected features to known true features and calculates 
TP, FP, FN, TN, sensitivity, and specificity.
}
\details{
Sensitivity = TP / (TP + FN). Specificity = TN / (TN + FP).
Use only features actually tested in `all` to get meaningful TNs.
}
\examples{
true <- c("gene1", "gene2", "gene3")
selected <- c("gene2", "gene3", "gene4")
all <- paste0("gene", 1:10)
compute_sens_spec(true, selected, all)

}

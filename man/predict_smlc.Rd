% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_predict_smlc.R
\name{predict_smlc}
\alias{predict_smlc}
\alias{predict_mlogit}
\title{Prediction based on the hidden genome sparse multinomial
logistic classifier}
\usage{
predict_smlc(fit, Xnew, Ynew = NULL, return_lin_pred = FALSE, ...)

predict_mlogit(fit, Xnew, Ynew = NULL, return_lin_pred = FALSE, ...)
}
\arguments{
\item{fit}{fitted hidden genome mlogit classifier, an output of fit_smlc.}

\item{Xnew}{test data design (or meta-design) matrix (observations
across rows and variables predictors/features across columns)
for which predictions are to be made from a fitted model. For a typical hidden
genome classifier this will be a matrix whose rows correspond to the test set
tumors, and columns correspond to (normalized by some functions of
the total mutation burdens in tumors) binary 1-0 presence/absence of
raw variants, counts of mutations at specific genes and counts of mutations
corresponding to specific mutation signatures etc.}

\item{Ynew}{the actual cancer categories for the test samples.
This is not used in computation, but is return as a component in the output,
for possibly easier post-processing.}
}
\value{
a list with entries (a) probs_predicted:
a \code{ncol(Xnew)} by n_cancer  (determined from \code{fit})
matrix of multinomial probabilities, providing
the predicted probability of each sample unit in Xnew
being classified into each cancer site,
and (b) predicted : a character vector listing
hard classes based on the predicted multinomial
probabilities (obtained by assigning individuals to
the classes with the highest predicted probabilities), and
optionally, (c) observed: if Ynew is supplied, then it
is returned as is.
}
\description{
Prediction based on the hidden genome sparse multinomial
logistic classifier
}
\note{
Predictors in \code{Xnew} that are not present in the
training set design matrix (stored in \code{fit}) are dropped, and predictors
not included in \code{Xnew} but present in training set design matrix are
all assumed to have zero values. This is convenient for a typical
hidden genome classifier where most predictors are (some normalized versions
of) counts (e.g. for gene and mutation signatures) or
binary presence/absence indicators (e.g., for raw variants) so that a zero
predictor value essentially indicates some form of "absence".
However, care must be taken for predictors whose 0 values
do not indicate absence.
}
\seealso{
fit_mlogit
}

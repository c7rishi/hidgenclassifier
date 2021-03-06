% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimal_threshold.R
\name{calc_one_v_rest_auc}
\alias{calc_one_v_rest_auc}
\title{Calculating area under Precision-Recall curve (PRC) and
Receiver-Operator characteristic curve (ROC) for all one-vs-rest
comparisons in the fitted model}
\usage{
calc_one_v_rest_auc(
  fit,
  measure = c("PRC", "ROC"),
  fitted_prob = NULL,
  include_baseline = TRUE,
  ...
)
}
\arguments{
\item{fit}{fitted hidden genome classifier object}

\item{measure}{Type of curve to use. Options include "PRC" (Precision Recall Curve) and
"ROC" (Receiver Operator characteristic Curve). Can be a vector.}

\item{fitted_prob}{an n_tumor x n_cancer matrix of predicted classification probabilities of
the \emph{training set tumors} to use for calculating ROC/PRC AUCs,
where n_tumor denotes the number of tumor/sample units,
and n_cancer is the number of cancer sites in the fitted hidden genome model (supplied
through \code{"fit"}). Row names and column names must
be identical to the the tumor/sample names and cancer labels used in  the fitted model. If \code{NULL}
(default) then the fitted probabilities are obtained from the model itself by either extracting pre-validated
predictive probabilities (only available for mlogit models), or simply using the fitted model to
make predictions on the training set.}

\item{include_baseline}{logical. Along with the computed \emph{observed} value(s) of the measure(s)
should the null baseline value(s) be returned. Here null baseline  refers to the expected
value of the corresponding measure associated with a "baseline" classifier that (uniform) randomly assigns
class labels to the sample units.}
}
\value{
Returns a data.table with \code{length(measure) + 1} columns
("Class" and measure(s)) and n_class + 1 many rows, where n_class
denotes the number of cancer types present in the fitted model; the
final row provides the Macro (average) metrics.
}
\description{
Calculating area under Precision-Recall curve (PRC) and
Receiver-Operator characteristic curve (ROC) for all one-vs-rest
comparisons in the fitted model
}
\details{
Under the hood, the function uses several functions from R package \code{precrec}
to compute the performance
metrics. The argument \code{fitted_prob}, when supplied, should ideally
contain predictive probabilities for training set tumors evaluated under a
cross-validation framework. If not supplied, pre-validated
prediction probabilities extracted from  mlogit models, and
overoptimistic prediction probabilities (obtained by simply using the fitted
model on the training data) for other models are used.
}

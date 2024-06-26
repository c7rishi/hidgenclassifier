% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_cancer_response.R
\name{extract_cancer_response}
\alias{extract_cancer_response}
\title{Convenience function for extracting cancer (response)
categories corresponding to all tumors
in a maf file for use in a hidden genome classifier}
\usage{
extract_cancer_response(
  maf,
  cancer_col = "cancer",
  sample_id_col = "sample",
  ...
)
}
\arguments{
\item{maf}{mutation annotation file --
a data frame-like object with at least two columns -- one providing
sample ids of tumor and one providing the associated cancer categories}

\item{cancer_col}{name of the column in \code{maf} that corresponds to cancer
sites for the tumor samples.}

\item{sample_id_col}{name of the column in \code{maf} containing tumor sample IDs.}

\item{...}{Unused.}
}
\value{
Returns a character vector containing cancer sites as determined
from cancer_col in maf, and named according to sample_id_col in maf.
}
\description{
Convenience function for extracting cancer (response)
categories corresponding to all tumors
in a maf file for use in a hidden genome classifier
}
\examples{
data("impact")
cancer_resp <- extract_cancer_response(
  maf = impact,
  cancer_col = "CANCER_SITE",
  sample_id_col = "patient_id"
)
head(cancer_resp)

}

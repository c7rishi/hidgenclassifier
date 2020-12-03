#' Convenience function for extracting cancer (response)
#' categories corresponding to all tumors
#' in a maf file for use in a hidden genome classifier
#'
#' @inheritParams variant_screen_mi
#'
#' @param maf mutation annotation file --
#' a data frame-like object with at least two columns -- one providing
#' sample ids of tumor and one providing the associated cancer categories
#'
#' @return
#'
#' Returns a character vector containing cancer sites as determined
#' from cancer_col in maf, and named according to sample_id_col in maf.
#'
#' @examples
#' data("impact")
#' cancer_resp <- extract_cancer_response(
#'   maf = impact,
#'   cancer_col = "CANCER_SITE",
#'   sample_id_col = "patient_id"
#' )
#' head(cancer_resp)
#'
#' @export
extract_cancer_response <- function(
  maf,
  cancer_col = "cancer",
  sample_id_col = "sample",
  ...
  ) {
  dt <- unique(
    data.table::as.data.table(
      maf
    )[,
      c(cancer_col, sample_id_col),
      with = FALSE
    ]
  )

  setNames(dt[[cancer_col]], dt[[sample_id_col]])
}

#' Extract Design-Metadesign (XU) matrix for a set of continuous
#' meta-features stored in a maf file
#' @inheritParams variant_screen_mi
#' @param maf mutation annotation file --
#' a data frame-like object containing columns for
#' variant labels, sample IDs and continuous meta-feature variables.
#' NOTE: uniqueness of rows of maf is assumed
#' @param mfeat_cols names of columns in \code{maf} containing continuous
#' meta-feature variables. Can be a vector.
#' @export

extract_design_mdesign_mcts <- function(
  maf,
  variant_col = "variant",
  mfeat_cols = "Histone_Marks",
  sample_id_col = "sample",
  ...
) {

  dt <- data.table::as.data.table(
    maf
  )[,
    c(sample_id_col, variant_col, mfeat_cols),
    with = FALSE
  ]

  setnames(
    dt,
    old = c(sample_id_col, variant_col),
    new = c("p", "v")
  )

  out <- dt[
    ,
    lapply(.SD, sum, na.rm = TRUE),
    .SDcols = mfeat_cols,
    by = p
  ] %>%
    magrittr::set_rownames(.$p) %>%
    .[, p := NULL] %>%
    data.matrix() %>%
    Matrix(sparse = TRUE)

  out
}

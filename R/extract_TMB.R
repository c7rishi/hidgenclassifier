#' #' Extract total mutation burden (TMB) of tumors in a maf file
#' @inheritParams extract_design
#'
#' @return
#' Returns a numeric vector providing the total
#' number of mutations observed in each tumor in maf, named
#' by the sample ids in maf.
#'
#' @examples
#' data("impact")
#' tmb <- extract_tmb(
#'   maf = impact,
#'   variant_col = "Variant",
#'   sample_id_col = "patient_id"
#' )
#' tmb
#'
#' @export
extract_tmb <- function(
  maf,
  variant_col = "variant",
  sample_id_col = "sample",
  ...
) {

  dt <- as.data.table(
    maf
  )[,
    c(sample_id_col, variant_col),
    with = FALSE
  ]

  setnames(
    dt,
    old = c(sample_id_col, variant_col),
    new = c("p", "v")
  )

  p <- v <- nv <- NULL # so that CRAN doesn't complain

  setkey(dt, p)

  dt_tmb <- dt[
    ,
    .(nv = length(unique(v))),
    by = p
  ]

  setNames(dt_tmb$nv, dt_tmb$p)
}

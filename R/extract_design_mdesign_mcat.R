#' Extract Design-Metadesign (XU) matrix for a categorical
#' meta-feature from a maf file
#' @inheritParams variant_screen_mi
#' @param maf mutation annotation file --
#' a data frame-like object with at least three columns containing
#' variant labels, sample IDs and (categorical) meta-feature labels.
#' NOTE: uniqueness of rows of maf is assumed.
#' @param mfeat_col name of the column in \code{maf} containing categorical
#' meta-feature labels.
#' @export
extract_design_mdesign_mcat <- function(
  maf,
  variant_col = "variant",
  mfeat_col = "gene",
  sample_id_col = "sample",
  mfeat_subset = NULL,
  ...
) {

  dt <- data.table::as.data.table(
    maf
  )[,
    c(sample_id_col, variant_col, mfeat_col),
    with = FALSE
  ]

  setnames(
    dt,
    old = c(sample_id_col, variant_col, mfeat_col),
    new = c("p", "v", "g")
  )

  data.table::setkey(dt, p, g)

  if (is.null(mfeat_subset)) {
    mfeat_subset <- unique(dt$g)
  }

  g_p <- unique(
    dt[
      g %in% mfeat_subset,
      .(nv = length(unique(v))),
      by = .(p, g)
    ]
  )[,
    `:=`(
      p1 = factor(p),
      g1 = factor(g)
    )
  ][,
    .(nv, p, p1, g, g1)
  ]


  Xmat <- Matrix::sparseMatrix(
    i = as.integer(g_p$p1),
    j = as.integer(g_p$g1),
    x = g_p$nv,
    dims = c(nlevels(g_p$p1), nlevels(g_p$g1)),
    dimnames = list(levels(g_p$p1), levels(g_p$g1))
  )

  all_p <- unique(dt$p)

  out <- Xmat %>%
    fill_sparsemat_zero(all_p, unique(mfeat_subset))

  out

}

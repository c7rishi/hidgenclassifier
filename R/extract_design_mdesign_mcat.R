#' Extract Design-Metadesign (XU) matrix for a categorical
#' meta-feature from a maf file
#' @inheritParams variant_screen_mi
#' @param maf mutation annotation file --
#' a data frame-like object with at least three columns containing
#' variant labels, sample IDs and (categorical) meta-feature labels.
#' NOTE: uniqueness of rows of maf is assumed.
#' @param mfeat_col name of the column in \code{maf} containing categorical
#' meta-feature labels.
#' @param mfeat_subset character vector providing the subset of categories of
#' the meta-feature for which the design matrix is to be created. If NULL
#' (default), all unique categories present in the \code{mfeat}_col in
#' \code{maf} is considered.
#'
#' @return
#' An n_tumor x n_gene sparse dgCMatrix, with (i, j)th entry providing the total
#' number of variants in tumor i associated with j-th meta-feature category,
#' as determined by \code{mfeat_col} of \code{maf}.
#'
#' @examples
#' data("impact")
#' gene_mdesign <- extract_design_mdesign_mcat(
#'   maf = impact,
#'   variant_col = "Variant",
#'   mfeat_col = "Hugo_Symbol",
#'   sample_id_col = "patient_id"
#' )
#' dim(gene_mdesign)
#'
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

  p <- v <- g <- NULL # so that CRAN check doesn't complain

  data.table::setkey(dt, p, g)

  if (is.null(mfeat_subset)) {
    mfeat_subset <- unique(dt$g)
  }

  nv <- p1 <- g1 <- NULL

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

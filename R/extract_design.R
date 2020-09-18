#' Obtain design matrix (X) for variants from a mutation annotation file (maf)
#'
#' @inheritParams variant_screen_mi
#' @param variant_subset The subset of variant for which the design matrix
#' is to be constructed. Either NULL, in which case all variants present in the
#' maf file will be considered, or a character vector which will constitute the
#' column names of the output.
#' @param maf mutation annotation file --
#' a data frame-like object with at least two columns containing
#' variant IDs and sample IDs. NOTE: uniqueness of rows of maf is assumed.
#' @return returns an n_tumor x n_variants design matrix for
#' variant indicators (in sparseMatrix format)
#'
#'
#' OUTPUT:
#' an N_patients x N_variants (or N_genes) (sparse) design
#' (or design-meta-design) matrix with appropriate row and
#' variant_col names.
#' @export
extract_design <- function(
  maf,
  variant_col = "variant",
  sample_id_col = "sample",
  variant_subset = NULL,
  ...
) {

  dt <- as.data.table(
    maf
  )[,
    c(variant_col, sample_id_col),
    with = FALSE
    ]

  setnames(
    dt,
    old = c(sample_id_col, variant_col),
    new = c("p", "vv")
  )

  setkey(dt, p, vv)

  if (is.null(variant_subset)) {
    vv_p <- unique(
      dt[,
         .(p, vv)
      ]
    )[,
      `:=`(
        p1 = factor(p),
        vv1 = factor(vv)
      )
    ][,
      .(p, p1, vv, vv1)
    ]

    variant_subset <- unique(dt$vv)
  } else {
    vv_p <- unique(
      dt[vv %in% variant_subset,
         .(p, vv)
      ]
    )[,
      `:=`(
        p1 = factor(p),
        vv1 = factor(vv)
      )
    ][,
      .(p, p1, vv, vv1)
    ]
  }

  Xmat <- Matrix::sparseMatrix(
    i = as.integer(vv_p$p1),
    j = as.integer(vv_p$vv1),
    x = 1,
    dims = c(nlevels(vv_p$p1), nlevels(vv_p$vv1)),
    dimnames = list(levels(vv_p$p1), levels(vv_p$vv1))
  )

  all_p <- unique(dt$p)

  out <- Xmat %>%
    fill_sparsemat_zero(all_p, unique(variant_subset))

  out
}

#' Mutual Information based feature screening of variants from a mutation
#' annotation file
#'
#'
#'
#' @param maf mutation annotation file --
#' a data frame-like object with at least three columns containing
#' variant labels, sample IDs, and cancer sites associated with the sample IDs.
#' NOTE: uniqueness of rows of maf is assumed.
#' @param variant_col name of the column in \code{maf} containing variant labels.
#' @param sample_id_col name of the column in \code{maf} containing tumor sample IDs.
#' @param cancer_col name of the column in \code{maf} that corresponds to cancer
#' sites for the tumor samples.
#' @param mi_rank_thresh rank threshold for screening variants. The top variants
#' with rank(MI_values) <= mi_rank_thresh is returned. Defaults to 250.
#' @param equal_cancer_prob_mi logical. Should the marginal probabilities of
#' cancer sites be assumed equal (i.e., uniform) while computing mutual
#' information? If \code{FALSE}, the relative frequencies of cancer sites in
#' maf are used. CAUTION: the (sample) relative frequencies of cancer sites
#' in \code{maf} may not necessarily be good approximations of the truth.
#' @param normalize_mi logical. Should mutual information be normalized by
#' product of square-roots of marginal Shannon entropies? Defaults to FALSE.
#' @param return_prob_mi logical. Should the computed mutual information and the
#' cancer site specific probabilities for these screened variants be returned?
#' Defaults to TRUE.
#' @param do_freq_screen logical. Should an overall (relative) frequency-based
#' screening be performed prior to MI based screening?
#' This may reduce the computation load substantially for whole genome
#' data where potentially tens of  millions of variants are observed only once. Defaults
#' to FALSE.
#' @param thresh_freq_screen Threshold for overall pan-cancer relative frequency
#' to use if a frequency-based screening is performed before mi based
#' screening. Defaults to 1/n_sample where n_sample is the pan-cancer
#' total number of tumors. Ignored if \code{do_freq_screen = FALSE}.
#'
#'
#' @param ... Unused.
#'
#' @details
#' The function first estimates via relative frequencies the cancer site
#' specific probabilities of encountering EACH variant in the maf file. Then using
#' these estimated probabilities and the marginal probabilities of cancer sites,
#' the (possibly normalized) mutual information between (a) the occurrence of a
#' variant-"j" in randomly chosen tumor and (b) the cancer site of the associated
#' tumor is computed for each variant-j in \code{maf}.
#' These MIs are then ranked and the variant labels associated with with
#' mi rank <= \code{mi_rank_thresh}  are returned.
#'
#' @return
#' a character vector listing the screened variant labels (sorted with the first
#' one having the highest MI) with ranks <= \code{mi_rank_thresh}.
#' Optionally, if \code{return_prob_mi = TRUE}, then
#' a data table named \code{prob_mi} listing cancer site specific probabilities
#' of ALL variants and the associated MIs are returned.
#'
#' @examples
#' data("impact")
#' top_v <- screen_variant_mi(
#'   maf = impact,
#'   variant_col = "Variant",
#'   cancer_col = "CANCER_SITE",
#'   sample_id_col = "patient_id",
#'   mi_rank_thresh = 200,
#'   return_prob_mi = FALSE
#' )
#' top_v
#'
#'
#' @export
screen_variant_mi <- function(
  maf,
  variant_col = "variant",
  cancer_col = "cancer",
  sample_id_col = "sample",
  equal_cancer_prob_mi = TRUE,
  return_prob_mi = TRUE,
  mi_rank_thresh = 250,
  normalize_mi = FALSE,
  do_freq_screen = FALSE,
  thresh_freq_screen = 1/length(unique(maf[[sample_id_col]])),
  ...
) {

  dots <- list(...)

  if (!is.null(dots$design_mat) &
      !is.null(dots$cancer_cat)) {


    # EXPERIMENTAL, uses a pre-computed design matrix

    Xmat <- Matrix(dots$design_mat, sparse = TRUE)

    dt_p_c <- data.table(
      p = names(dots$cancer_cat),
      c = unname(dots$cancer_cat)
    )

    # number of tumors per cancer site
    np_c <- dt_p_c[
      ,
      .(np = length(unique(p))),
      by = c
    ][,
      prop := ifelse(
        equal_cancer_prob_mi,
        1/.N,
        np/sum(np)
      )
    ]

    v_rf_c_mi <- dt_p_c[
      ,
      {
        as.list(Matrix::colMeans(Xmat[.SD$p, ]))
      },
      by = c
    ] %>%
      melt(
        id.vars = "c",
        variable.name = "vv",
        value.name = "v_rf"
      ) %>%
      dcast(
        vv ~ c,
        fill = 0,
        value.var = "v_rf"
      ) %>%
      .[, vv := as.character(vv)]


  } else {

    dt <- as.data.table(
      maf
    )[,
      c(variant_col, cancer_col, sample_id_col),
      with = FALSE
    ]

    setnames(
      dt,
      old = c(variant_col, cancer_col, sample_id_col),
      new = c("vv", "c", "p")
    )
    setkey(dt, vv, c)

    # number of tumors per cancer site
    np_c <- dt[
      ,
      .(np = length(unique(p))),
      by = c
    ][,
      prop := ifelse(
        equal_cancer_prob_mi,
        1/.N,
        np/sum(np)
      )
    ]

    # add cancer specific-sample size column
    # per cancer site to the maf
    dt[,
       np := length(unique(p)),
       by = c
    ]

    if (do_freq_screen) {
      np_full <- sum(np_c$np)
      v_rf_marg <- dt[
        ,
        .(v_rf = length(unique(p))/np_full),
        by = vv
      ]

      v_keep <- v_rf_marg[v_rf > thresh_freq_screen]$vv

      dt <- dt[vv %in% v_keep]
    }


    v_rf_c_mi <- dt[
      # compute relative frequencies of
      # variants by cancer
      # then dcast
      ,
      {
        freq <- length(unique(p))
        .(v_rf = freq/np[1],
          v_f = freq)
      },
      by = .(vv, c)
    ][,
      is_singleton := (max(v_f) == 1),
      by = vv
    ][
      is_singleton == FALSE
    ] %>%
      dcast(
        vv ~ c, fill = 0,
        value.var = "v_rf"
      )

  }

  cancer_prob <- np_c$prop
  names(cancer_prob) <- np_c$c

  # add mi values for the variants
  # along with its (descending) ranks
  v_rf_c_mi[
    ,
    mi := calc_minfo(
      data.matrix(.SD),
      cancer_prob,
      normalize = normalize_mi
    ),
    .SDcols = names(cancer_prob)
  ][,
    mi_rank := rank(-mi)
  ]


  # arrange by descending mi
  setorder(v_rf_c_mi, -mi)
  setnames(v_rf_c_mi, "vv", variant_col)


  out_mat <- v_rf_c_mi[mi_rank <= mi_rank_thresh]
  out <- as.character(out_mat[[variant_col]])

  if (return_prob_mi) {
    attr(out, "prob_mi") <- out_mat
  }

  out
}


#' @rdname screen_variant_mi
#' @export
variant_screen_mi <- screen_variant_mi

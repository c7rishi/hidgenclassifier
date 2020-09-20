#' Find design (or design-meta-design) matrix for
#' variants and genes
#'
#' @details
#' INPUT:
#' maf: mutation annotation file (maf) with colum names
#' g = gene, v = variants, c = cancer type,
#' and p = patient (tumor) id.
#' column: name of the column that corresponds
#' to variants. Either "v" (for variants, if design matrix) or
#' "g" for genes (if design-meta-design)
#' col_select: either "all", which means all columns of the design
#' matrix will be retained, based on all variants/genes observed in maf,
#' or a character vector of variants/genes that are to be retained.
#'
#' OUTPUT:
#' an N_patients x N_variants (or N_genes) (sparse) design
#' (or design-meta-design) matrix with appropriate row and
#' column names.
#' @export
make_design <- function(maf,
                        column = "v",
                        sampleid_col = "p",
                        col_select = "all",
                        calc_nontarget_p = TRUE,
                        target_g = "all") {

  fill_sparsemat_zero <- function(mat, rownames, colnames) {
    # mat_sparse <- Matrix::Matrix(mat, sparse = TRUE)
    zeromat <- Matrix::Matrix(data = 0,
                              nrow = length(rownames),
                              ncol = length(colnames),
                              dimnames = list(rownames, colnames),
                              sparse = TRUE)
    rownames_common <- intersect(rownames, rownames(mat))
    colnames_common <- intersect(colnames, colnames(mat))
    zeromat[rownames_common, colnames_common] <-
      mat[rownames_common, colnames_common]
    zeromat[rownames, colnames]
  }


  `%>%` <- magrittr::`%>%`
  dt <- data.table::as.data.table(maf)[, vv := get(column)]
  # browser()

  data.table::setnames(
    dt,
    old = c(sampleid_col),
    new = c("p")
  )

  data.table::setkey(dt, p, vv)


  if (all(col_select == "all")) {
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

    col_select <- unique(dt[[column]])
  } else {
    vv_p <- unique(
      dt[vv %in% col_select,
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
  if (nrow(Xmat) < length(all_p)) {
    Xmat <- fill_sparsemat_zero(Xmat, all_p, colnames(Xmat))
  }


  if (all(target_g == "all")) {
    g_nontarget <- NULL
  } else {
    g_nontarget <- unique(dt$g) %>%
      setdiff(target_g)

    if (length(g_nontarget) > 0) {
      nsamp <- length(unique(dt$p))
      pv <- dt[
        g %in% g_nontarget &
          v %in% col_select
      ][,
        .(v_rf = length(unique(p))/nsamp),
        by = v
      ]

      p.x <- pv$v_rf
      names(p.x) <- pv$v

      attr(Xmat, "p_v") <- p.x
    }

  }


  Xmat
}


#' @export
extract_design_metadesign_metacat <- function(maf,
                                              metafeat_col = "g",
                                              metafeat_levels = "all",
                                              calc_nontarget_p.U = TRUE,
                                              target_g = "all",
                                              p_list = "all"
) {

  `%>%` <- magrittr::`%>%`
  dt <- data.table::as.data.table(maf)
  data.table::setkey(dt, p)

  dt$metacol <- dt[[metafeat_col]]

  if (all(target_g == "all")) {
    target_g <- unique(dt$g)
  }

  if (all(metafeat_levels == "all")) {
    metafeat_levels <- unique(dt$metacol)
  }

  if (all(p_list == "all")) {
    p_list <- unique(dt$p)
  }

  metafeat_levels <- intersect(dt$metacol, metafeat_levels)

  data.table::setkey(dt, p, metacol)
  p_p1 <- unique(dt[, .(p)])[, p1 := 1:.N]
  meta_meta1 <- unique(
    dt[g %in% target_g &
         metacol %in% metafeat_levels,
       .(metacol)
    ]
  )[, metacol1 := 1:.N]

  nv_by_meta_p <- merge(dt[g %in% target_g, .(p, metacol, v)], p_p1, by = "p") %>%
    merge(meta_meta1, by = "metacol") %>%
    .[metacol %in% metafeat_levels,
      .(
        n_v = length(unique(v))
      ),
      by = .(p1, metacol1)
    ]

  XUmat <- Matrix::sparseMatrix(
    i = nv_by_meta_p$p1,
    j = nv_by_meta_p$metacol1,
    x = nv_by_meta_p$n_v,
    dims = c(max(p_p1$p1), max(meta_meta1$metacol1))
  )

  colnames(XUmat) <- meta_meta1$metacol
  rownames(XUmat) <- p_p1$p

  XUmat_adj <- Matrix::Matrix(
    0,
    nrow = length(p_list),
    ncol = length(metafeat_levels),
    dimnames = list(p_list, metafeat_levels)
  )
  XUmat_adj[rownames(XUmat), colnames(XUmat)] <- XUmat

  #
  # g_non_target <- unique(dt$g) %>%
  #   setdiff(target_g)
  #
  #   # interset with selected (important) columns,
  #   # if metafeat_col = "g
  #   if (metafeat_col == "g") {
  #     g_non_target <- g_non_target %>%
  #       intersect(metafeat_levels)
  #   }
  #
  #
  # if (length(g_non_target) > 0) {
  #   nsamp <- length(unique(dt$p))
  #   pv_metacol <- dt[
  #     g %in% g_non_target
  #     ][,
  #       .(
  #         metacol = metacol[1],
  #         v_rf = length(unique(p))/nsamp
  #       ),
  #       by = v
  #       ][,
  #         .(pv_sum = sum(v_rf)),
  #         by = metacol
  #         ]
  #
  #   nontarget_p.U <- pv_metacol$pv_sum
  #   names(nontarget_p.U) <- pv_metacol$metacol
  #
  #   # p.U <- nontarget_p.U[metafeat_levels] %>%
  #   #   ifelse(is.na(.), 0, .)
  #
  #   attr(XUmat, "p.U") <- nontarget_p.U
  # }

  XUmat_adj
}

#' @export
extract_design_metadesign_cts <- function(maf,
                                          meta_cols = cols) {
  dt <- data.table::as.data.table(maf)
  data.table::setkey(dt, p)

  XU_dt <- dt[,
              lapply(.SD, sum, na.rm = TRUE),
              keyby = p,
              .SDcols = meta_cols
  ]

  XUmat <- data.matrix(XU_dt[, meta_cols, with = FALSE])
  rownames(XUmat) <- XU_dt$p

  XUmat
}

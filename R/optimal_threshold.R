#' Determining optimal one-vs-rest classification
#' thresholds from fitted hidden genome models using
#' predicting performance measures
#'
#' @param fit fitted hidden genome classifier object
#' @param measure prediction assessment measure.
#'
#' @return
#' If \code{length(measure) == 1} the function returns a named vector with optimal
#' one-vs-rest classification thresholds for all cancer
#' classes in the fitted hidden genome object (fit). The optimal
#' values obtained at the corresponding optimal thresholds are
#' returned as an attribute "optimal_value".
#'
#' If \code{length(measure) > 1} a named list is returned, with each
#' entry providing the optimal thresholds across all cancer categories
#' (along with the associated optimal measure values as an attribute)
#'corresponding to that \code{measure}.
#' @export

optimal_threshold <- function(fit,
                              measure = "fscore",
                              ...) {

  meth <- fit$method



  # if (any(!measure %in% ))

  if (!fit$method == "mlogit") {
    stop("Only implemented for mlogit classifier")
  }

  if (fit$method == "mlogit") {

    probs_predicted <- fit$fit$fit.preval[, , "s0"] %>%
      apply(1, softmax) %>%
      t()
  }

  resp_class <- colnames(probs_predicted)

  probs_pred_df <- probs_predicted %>%
    data.matrix() %>%
    data.table::data.table(keep.rownames=TRUE) %>%
    data.table::setnames("rn", "pid") %>%
    .[,
      `:=`(
        obs_canc = fit$Y,
        dsid = 1,
        meth = meth
      )
    ] %>%
    data.table::melt(
      measure.vars = resp_class,
      value.name = "pred_prob",
      variable.name = "canc"
    ) %>%
    .[, obs_indic := as.numeric(obs_canc == canc)] %>%
    # .[, obs_canc := NULL] %>%
    unique()

  all_classes <- colnames(probs_predicted)

  thresh_res <- lapply(
    all_classes,
    function(this_class) {
      indiv_one_v_rest_comparison(
        probs_pred_df[canc == this_class],
        measure = measure
      )[,
        class := this_class
      ][,
        .(class, measure,
          optimal_threshold, optimal_value)
      ]
    }
  ) %>%
    do.call(rbind, .)

  out <- sapply(
    measure,
    function(this_meas) {
      tmp <- thresh_res[measure == this_meas]

      out_thresh <- tmp$optimal_threshold %>%
        setNames(tmp$class)
      out_val <- tmp$optimal_value %>%
        setNames(tmp$class)
      attr(out_thresh, "optimal_value") <- out_val

      out_thresh
    },
    USE.NAMES = TRUE,
    simplify = FALSE
  )

  if (length(measure) == 1) {
    out <- out[[1]]
  }

  out
}




indiv_one_v_rest_comparison <- function(
  df,
  measure = c("fscore", "mcc"),
  ...) {
  unique_dsids <- unique(df$dsid) %>% sort()

  # df <- unique(df)

  scores_list <- lapply(
    unique_dsids,
    function(this_dsid) {

      df[dsid == this_dsid,
         .(pid, meth, pred_prob)] %>%
        data.table::dcast(
          pid ~ meth,
          value.var = "pred_prob"
        ) %>%
        data.table::setorder(pid) %>%
        .[, pid := NULL] %>%
        data.table::setDF() %>%
        as.list()
    }
  )
  scores_list_len <- sapply(scores_list, length)
  keep_dsid <- scores_list_len == max(scores_list_len)
  # keep_meth <-
  scores_list <- scores_list[keep_dsid]
  unique_dsids <- unique_dsids[keep_dsid]

  n_meth <- length(scores_list[[1]])
  all_meth <- names(scores_list[[1]])

  labels <- data.table::copy(df)[
    ,
    ndsid_meth := length(unique(dsid)),
    by = meth
  ][
    dsid == min(dsid[keep_dsid])
  ][
    ndsid_meth == max(ndsid_meth)
  ][
    meth == tail(unique(meth), 1)
  ] %>%
    data.table::setorder(pid) %>%
    .$obs_indic

  dsids <- unique_dsids %>%
    lapply(
      function(this_dsid) {
        rep(this_dsid, n_meth)
      }
    ) %>%
    unlist()

  modnames <- scores_list %>%
    lapply(names) %>%
    unlist()


  mmdata_obj <- precrec::mmdata(
    scores = scores_list,
    labels = labels,
    dsids = dsids,
    posclass = 1,
    modnames = modnames
  )

  performance <- mmdata_obj %>%
    precrec::evalmod(mode = "basic")

  out <- lapply(
    measure,
    function(this_measure) {
      get_thresh_mmevalmod(
        mmdata_obj = mmdata_obj,
        evalmod_obj = performance,
        measure_name = this_measure
      ) %>%
        .[,
          `:=`(
            measure = this_measure,
            modname = NULL
          )
        ] %>%
        data.table::setnames(
          c(this_measure, "opt_thresh"),
          c("optimal_value", "optimal_threshold")
        )
    }
  ) %>%
    do.call(rbind, .)

  out
}



get_thresh_mmevalmod <- function(mmdata_obj,
                                 evalmod_obj,
                                 measure_name = "fscore") {

  out <- autoplot(
    evalmod_obj,
    measure_name,
    show_cb = FALSE,
    plot = FALSE
  )$data %>%
    data.table::as.data.table() %>%
    # grouped filter
    .[, .SD[y == max(y, na.rm = TRUE)[1]], by = modname] %>%
    .[, .(modname, x, y)] %>%
    data.table::setnames(
      c("x", "y"),
      c("normalized_rank", measure_name)
    )

  # associate threshold with normalized ranks
  score_rank_list <- mmdata_obj %>%
    lapply(
      # precrec::evalmod(mode = "basic") %>%
      function(yy) {
        yy[c("scores", "ranks")] %>%
          as.list() %>%
          data.frame() %>%
          data.table::setDT() %>%
          .[,
            normalized_rank := (ranks - 1)/(length(ranks) -1)
          ] %>%
          data.table::setorder(normalized_rank) %>%
          .[, .(normalized_rank, scores)] %>%
          unique()
      }
    ) %>%
    setNames(
      sapply(mmdata_obj, attr, "modname")
    )

  score_rank_fn_list <- lapply(
    score_rank_list,
    function(this_score_rank) {
      approxfun(
        x = this_score_rank$normalized_rank,
        y = this_score_rank$scores
      )
    }
  )

  out$opt_thresh <- mapply(
    function(this_normalized_rank, this_score_rank_fun) {
      this_score_rank_fun(this_normalized_rank)
    },
    out[[measure_name]],
    score_rank_fn_list[out$modname],
    SIMPLIFY = FALSE
  ) %>%
    setNames(out$modname)

  out[, normalized_rank := NULL]

  out
}

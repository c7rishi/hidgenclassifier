# TODO: add an argument "Y" or "true_labels" to
# calc_one_v_rest_auc so that these can be computed
# without requiring a single fitted hidden genome model
# as input (e.g., when predictions are obtained from
# multiple models)


#' Determine optimal one-vs-rest classification
#' thresholds from fitted hidden genome models using
#' prediction-based performance measures
#'
#' @param fit fitted hidden genome classifier object
#' @param measure prediction assessment measure. Options include "fscore",
#' "mcc" (Mathews Correlation Coefficient). Can be a vector.
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
#'
#' @export
optimal_threshold <- function(fit,
                              measure = "fscore",
                              fitted_prob = NULL,
                              true_labels = NULL,
                              ...) {

  meth <- fit$method

  stopifnot(requireNamespace("precrec"))

  probs_pred_df <- create_pred_prob_df_from_fit(fit, fitted_prob)

  if (attr(probs_pred_df, "msg") != "") {
    warning(attr(probs_pred_df, "msg"))
  }

  all_classes <- names(fit$alpha)

  class <- NULL # so that R check doesn't complain

  thresh_res <- lapply(
    all_classes,
    function(this_class) {
      indiv_one_v_rest_hard_comparison(
        probs_pred_df[canc == this_class],
        measure = measure
      )[,
        class := this_class
      ][,
        list(class, measure,
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



#' Calculating area under Precision-Recall curve (PRC) and
#' Receiver-Operator characteristic curve (ROC) for all one-vs-rest
#' comparisons in the fitted model
#' @inheritParams optimal_threshold
#' @param measure Type of curve to use. Options include "PRC" (Precision Recall Curve) and
#' "ROC" (Receiver Operator characteristic Curve). Can be a vector.
#' @param fitted_prob  an n_tumor x n_cancer matrix of predicted classification probabilities of
#' the *training set tumors* to use for calculating ROC/PRC AUCs,
#' where n_tumor denotes the number of tumor/sample units,
#' and n_cancer is the number of cancer sites in the fitted hidden genome model (supplied
#' through \code{"fit"}). Row names and column names must
#' be identical to the the tumor/sample names and cancer labels used in  the fitted model. If \code{NULL}
#' (default) then the fitted probabilities are obtained from the model itself by either extracting pre-validated
#' predictive probabilities (only available for mlogit models), or simply using the fitted model to
#' make predictions on the training set.
#' @param include_baseline logical. Along with the computed *observed* value(s) of the measure(s)
#'  should the null baseline value(s) be returned. Here null baseline  refers to the expected
#'  value of the corresponding measure associated with a "baseline" classifier that (uniform) randomly assigns
#'  class labels to the sample units.
#'
#'
#' @details
#' Under the hood, the function uses several functions from R package \code{precrec}
#' to compute the performance
#' metrics. The argument \code{fitted_prob}, when supplied, should ideally
#' contain predictive probabilities for training set tumors evaluated under a
#' cross-validation framework. If not supplied, pre-validated
#' prediction probabilities extracted from  mlogit models, and
#' overoptimistic prediction probabilities (obtained by simply using the fitted
#' model on the training data) for other models are used.
#'
#'
#' @return
#' Returns a data.table with \code{length(measure) + 1} columns
#' ("Class" and measure(s)) (\code{2 * length(measure) + 1} many columns if
#' \code{include_baseline = TRUE}) and n_class + 1 many rows, where n_class
#' denotes the number of cancer types present in the fitted model; the
#' final row provides the Macro (average) metrics.
#' @export
calc_one_v_rest_auc <- function(fit,
                                measure = c("PRC", "ROC"),
                                fitted_prob = NULL,
                                include_baseline = TRUE,
                                ...) {

  meth <- fit$method

  if (any(!measure %in% c("PRC", "ROC"))) {
    stop('"measure" must either be "PRC" or "ROC" or both')
  }

  stopifnot(requireNamespace("precrec"))

  probs_pred_df <- create_pred_prob_df_from_fit(fit, fitted_prob)

  if (attr(probs_pred_df, "msg") != "") {
    warning(attr(probs_pred_df, "msg"))
  }

  all_classes <- probs_pred_df$canc %>%
    as.character() %>%
    unique() %>%
    sort()

  class <- PRC <- ROC <- NULL # so that R check doesn't complain

  all_measure <- c("PRC", "ROC")

  out <- lapply(
    all_classes,
    function(this_class) {
      indiv_one_v_rest_soft_comparison(
        probs_pred_df[canc == this_class],
        measure = all_measure
      )[,
        class := this_class
      ][,
        c("class", toupper(all_measure)),
        with = FALSE
      ]
    }
  ) %>%
    do.call(rbind, .) %>%
    rbind(
      .[,
        list(
          class = "Macro (Average)",
          PRC = mean(PRC, na.rm = TRUE),
          ROC = mean(ROC)
        )
      ]
    )

  prc_baseline <- fit$Y %>%
    unname() %>%
    table() %>%
    c() %>%
    {./sum(.)} %>%
    c("Macro (Average)" = mean(.))

  if (include_baseline) {
    out[
      ,
      `:=`(
        PRC_baseline = prc_baseline[class],
        ROC_baseline = 0.5
      )
    ]
  }

  final_colnames <- names(out) %>%
    grep(
      paste(measure, collapse = "|"), .,
      ignore.case = TRUE, value = TRUE
    ) %>%
    {c("class", .)}


  out[, final_colnames, with = FALSE]

}

# processed prediction data table (to use in precrec)
# using prediction probability matrix
create_pred_prob_df_from_fit <- function(fit, fitted_prob) {

  meth <- fit$method

  if (is.null(fitted_prob)) {
    probs_predicted_mat <- create_pred_prob_matrix_from_fit(fit)
  } else {
    probs_predicted_mat <- fitted_prob
    # TODO: insert checks for fitted probabilities
    attr(probs_predicted_mat, "msg") <- ""
  }

  resp_class <- colnames(probs_predicted_mat)

  obs_canc <- dsid <- NULL # so that R check doesn't complain

  probs_pred_df <- probs_predicted_mat %>%
    data.matrix() %>%
    data.table::data.table(keep.rownames = TRUE) %>%
    data.table::setnames("rn", "pid") %>%
    .[,
      `:=`(
        obs_canc = fit$Y[rownames(probs_predicted_mat)],
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

  attr(probs_pred_df, "msg") <- attr(probs_predicted_mat, "msg")

  probs_pred_df
}


# create predicted probability matrix
# (prevalidated for mlogit, and simple overoptimistic
# prediction for other models)
# for the training set individuals
create_pred_prob_matrix_from_fit <- function(fit) {
  meth <- fit$method

  if (meth == "mlogit" & !is.null(fit$fit$fit.preval)) {
    cvf <- fit$fit

    probs_predicted <- cvf$fit.preval[, , cvf$lambda == cvf$lambda.1se] %>%
      apply(1, softmax) %>%
      t()

    msg <- ""

  } else {
    Xnew <- fit$X
    pred_obj <- do.call(
      paste0("predict_", meth),
      list(fit = fit, Xnew = Xnew)
    )
    probs_predicted <- pred_obj$probs_predicted

    msg <- if (meth == "mlogit") {
      msg <- paste(
        "Prevalidated prediction probabilities for training samples",
        "are missing in the fitted mlogit model.",
        "Overoptimistic prediction probabilities for training samples are used.",
        "Either ensure 'keep = TRUE' (default) when running fit_mlogit",
        "or supply cross-validation prediction probabilities for the training samples",
        "via the argument 'fitted_prob'",
        collapse = " "
      )
    } else {
      msg <- paste(
        paste0(
          "Prevalidated prediction probabilities for training samples",
          "are not available for", meth, "models."
        ),
        "Overoptimistic prediction probabilities for training samples are used.",
        "upply cross-validation prediction probabilities for the training samples",
        "via the argument 'fitted_prob'",
        collapse = " "
      )
    }


  }

  attr(probs_predicted, "msg") <- msg

  probs_predicted
}


get_mmdata_obj <- function(df, ...) {

  stopifnot(requireNamespace("precrec"))

  unique_dsids <- unique(df$dsid) %>% sort()

  dsid <-pid <- meth <- pred_prob <- NULL

  scores_list <- lapply(
    unique_dsids,
    function(this_dsid) {

      df[dsid == this_dsid,
         list(pid, meth, pred_prob)] %>%
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

  ndsid_meth <- NULL
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

  mmdata_obj
}


indiv_one_v_rest_hard_comparison <- function(
  df,
  measure = c("fscore", "mcc"),
  ...) {


  stopifnot(requireNamespace("precrec"))

  mmdata_obj <- get_mmdata_obj(df)

  performance <- mmdata_obj %>%
    precrec::evalmod(mode = "basic")

  measure <- modname <- NULL

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

  .SD <- data.table::.SD

  modname <- x <- y <- NULL  # so that R check doesn't complain

  out <- autoplot(
    evalmod_obj,
    measure_name,
    show_cb = FALSE,
    plot = FALSE
  )$data %>%
    data.table::as.data.table() %>%
    # grouped filter
    .[, .SD[y == max(y, na.rm = TRUE)[1]], by = modname] %>%
    .[, list(modname, x, y)] %>%
    data.table::setnames(
      c("x", "y"),
      c("normalized_rank", measure_name)
    )

  normalized_rank <- ranks <- scores <- NULL

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
          .[, list(normalized_rank, scores)] %>%
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

indiv_one_v_rest_soft_comparison <- function(df,
                                             measure = c("PRC", "ROC")) {
  stopifnot(requireNamespace("precrec"))

  mmdata_obj <- get_mmdata_obj(df)

  modnames <- curvetypes <- NULL

  out <- mmdata_obj %>%
    precrec::evalmod(mode = "prcroc") %>%
    precrec::auc() %>%
    data.table::setDT() %>%
    data.table::dcast(
      modnames ~ curvetypes,
      value.var = "aucs"
    ) %>%
    .[,
      c("modnames", toupper(measure)),
      with = FALSE
    ]

  out
}

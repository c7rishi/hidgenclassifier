as.data.frame <- function(x) {
  if (is(class(x), "dgCMatrix")) {
    as.data.frame(data.matrix(x))
  } else {
    as.data.frame(x)
  }
}

#' Hidden genome random forest classifier (rfc)
#'
#' @details Light wrapper around randomForest or ranger to use in hidden
#' genome classification
#' @param backend Which backend to use? Available options are
#' "ranger" and "randomForest" corresponding to the respective R packages.
#' NOTE: randomForest does not support sparseMatrix, and the predictor matrix
#' is coerced into an ordinary matrix. This means using randomForest will likely
#' be more memory intensive and hence slower than ranger.
#' @param ... additional arguments passed to ranger::ranger or randomForest::randomForest (depending  on backend).
#' @param tune logical. Tune the random forest hyper parameters? Only used if
#' backend = "ranger". Defaults to TRUE. If TRUE, a list of models are trained with
#' various mtry and num.trees parameters, and the fitted model with minimum oob
#' prediction error is returned.
#' @inheritParams fit_smlc
#'
#' @examples
#' data("impact")
#' top_v <- variant_screen_mi(
#'   maf = impact,
#'   variant_col = "Variant",
#'   cancer_col = "CANCER_SITE",
#'   sample_id_col = "patient_id",
#'   mi_rank_thresh = 50,
#'   return_prob_mi = FALSE
#' )
#' var_design <- extract_design(
#'   maf = impact,
#'   variant_col = "Variant",
#'   sample_id_col = "patient_id",
#'   variant_subset = top_v
#' )
#'
#' canc_resp <- extract_cancer_response(
#'   maf = impact,
#'   cancer_col = "CANCER_SITE",
#'   sample_id_col = "patient_id"
#' )
#' pid <- names(canc_resp)
#' # create five stratified random folds
#' # based on the response cancer categories
#' set.seed(42)
#' folds <- data.table::data.table(
#'   resp = canc_resp
#' )[,
#'   foldid := sample(rep(1:5, length.out = .N)),
#'   by = resp
#' ]$foldid
#'
#' # 80%-20% stratified separation of training and
#' # test set tumors
#' idx_train <- pid[folds != 5]
#' idx_test <- pid[folds == 5]
#'
#' # train a classifier on the training set
#' # using only variants (will have low accuracy
#' # -- no meta-feature information used)
#' fit0 <- fit_rfc(
#'   X = var_design[idx_train, ],
#'   Y = canc_resp[idx_train],
#'   tune = FALSE
#' )
#'
#' pred0 <- predict_rfc(
#'   fit = fit0,
#'   Xnew = var_design[idx_test, ]
#' )
#'
#' @export
fit_rfc <- function(
  X, Y, backend = "ranger",
  tune = TRUE,
  mtry = NULL,
  n_mtry = 6,
  max.depth = c(0, 10^(-4:1)),
  num.trees = 1000,
  ...
) {
  if (!backend %in% c("randomForest", "ranger")) {
    stop('backend must be one of "randomForest" or "ranger"')
  }

  dots <- list(...)

  if (backend == "ranger") {
    if (!tune) {
      fit <- ranger::ranger(
        x = X,
        y = as.factor(Y),
        classification = TRUE,
        probability = TRUE,
        mtry = mtry,
        num.trees = num.trees,
        max.depth = max.depth[1],
        ...
      )
    } else {

      n_X <- ncol(X)
      if (is.null(mtry)) {
        mtry <- seq(
          floor(n_X^0.3),
          floor(n_X^0.7),
          length.out = n_mtry
        ) %>%
          floor()
      }

      inparam_list <- expand.grid(
        mtry = mtry,
        max.depth = max.depth
      )

      fit_list <- mapply(
        function(this_mtry, this_max.depth) {
          ranger::ranger(
            x = X,
            y = as.factor(Y),
            classification = TRUE,
            probability = TRUE,
            num.trees = num.trees,
            mtry = this_mtry,
            max.depth = this_max.depth,
            ...
          )
        },
        this_mtry = inparam_list$mtry,
        this_max.depth = inparam_list$max.depth,
        SIMPLIFY = FALSE
      )

      oob_error <- sapply(fit_list, "[[", "prediction.error")
      fit <- fit_list[[which.min(oob_error)[1]]]

    }
  }
  else {
    fit <- randomForest::randomForest(
      x = as.matrix(X),
      y = as.factor(Y),
      ...
    )
  }
  out <- list(
    X = X,
    Y = Y,
    fit = fit,
    backend = backend,
    method = "rf"
  )

  out
}




#' prediction based on hidden genome random forest classifier
#' @param fit Fitted random forest hidden genome classifier (output of
#' fit_rfc).
#' @export
predict_rfc <- function(fit,
                        Xnew,
                        Ynew = NULL, ...) {
  fit_rf <- fit$fit

  if (fit$backend == "ranger") {
    Xold_names <- fit_rf$forest$independent.variable.names
  } else {
    Xold_names <- rownames(fit_rf$importance)
  }

  # Xnew_adj <- adjust_Xnew(Xnew, Xold_names)
  Xnew_adj <- Xnew %>%
    fill_sparsemat_zero(
      rownames = rownames(.),
      colnames = Xold_names
    )

  if (fit$backend == "ranger") {
    predict_obj <- predict(
      fit_rf,
      data = Xnew_adj,
      type = "response"
    )

    predict_prob <- predict_obj$predictions %>%
      magrittr::set_rownames(rownames(Xnew_adj)) %>%
      .[, sort(colnames(.))]
  } else {
    Xnew_adj <- as.matrix(Xnew_adj)
    predict_prob <- as.matrix(
      predict(
        fit_rf,
        newdata = as.matrix(Xnew_adj),
        type = "prob"
      )
    )[rownames(Xnew_adj), ] %>%
      .[, sort(colnames(.))]
  }

  pred_class <- apply(
    predict_prob,
    1,
    function(x) names(x)[which.max(x)]
  )

  list(
    "predicted" = pred_class,
    "probs_predicted" = predict_prob,
    "observed" = Ynew
  )
}

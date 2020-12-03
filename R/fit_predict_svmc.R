#' Hidden genome SVM classifier (svmc)
#'
#' @details Light wrapper around e1071::svm or liquidSVM::mcSVM to use in hidden
#' genome classification
#' @param ... additional arguments passed to e1071:tune.svm, or
#' liquidSVM::svm.
#' @param backend the backend to use. Either "e1071" or "liquidSVM". Defaults to
#' "liquidSVM"
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
#' fit0 <- fit_svmc(
#'   X = var_design[idx_train, ],
#'   Y = canc_resp[idx_train]
#' )
#'
#' pred0 <- predict_svmc(
#'   fit = fit0,
#'   Xnew = var_design[idx_test, ]
#' )
#'
#'
#' @export
fit_svmc <- function(X,
                     Y,
                     backend = "liquidSVM",
                     scale = TRUE,
                     scale_fn = function(x) 2*sd(x),
                     ...) {
  dots <- list(...)
  if (is.null(dots$gamma)) {
    gamma <- 10^(-4:3)
  }
  if (is.null(dots$cost)) {
    dots$cost <- 10^(-4:3)
  }


  if (scale) {

    X_scale <- apply(
      X, 2,
      match.fun(scale_fn)
    ) %>%
      pmax(1) %>%
      setNames(colnames(X))

    X <- X %>%
      divide_cols(X_scale)

    attr(X, "scaled:scale") <- X_scale
  }

  if (backend == "e1071") {
    fit <- list(
      x = X,
      y = as.factor(Y),
      probability = TRUE,
      kernel = "radial"
    ) %>%
      c(dots) %>%
      do.call(e1071::tune.svm, .)
  } else if (backend == "liquidSVM") {

    if (is.null(dots$type)) {
      dots$type <- "AvA_ls"
    }

    if (is.null(dots$max_gamma)) {
      dots$max_gamma <- 1e6
    }

    if (is.null(dots$min_gamma)) {
      dots$min_gamma <- 1e-8
    }

    if (is.null(dots$min_lambda)) {
      dots$min_lambda <- 1e-8
    }

    fit <- list(x = as.matrix(X),
                y = as.factor(Y),
                predict.prob = TRUE) %>%
      c(dots) %>%
      do.call(liquidSVM::mcSVM, .)

  }
  out <- list(
    X = X,
    Y = Y,
    fit = fit,
    scale = scale,
    backend = backend,
    method = "svm"
  )

  out
}




#' prediction based on hidden genome random forest classifier
#' @seealso fit_svmc
#' @inheritParams predict_mlogit
#' @param fit fitted hidden genome SVM classifier (output of
#' \code{fit_svmc()})
#' @export
predict_svmc <- function(fit,
                         Xnew,
                         Ynew = NULL, ...)  {

  Xold_names <- colnames(fit$X)
  Xnew_adj <- Xnew %>%
    fill_sparsemat_zero(
      rownames = rownames(.),
      colnames = Xold_names
    )


  if (fit$scale) {
    Xscale <- rep(1, ncol(Xnew_adj)) %>%
      setNames(colnames(Xnew_adj))

    Xscale[Xold_names] <- attr(
      fit$X,
      "scaled:scale"
    )[Xold_names]

    Xnew_adj <- Xnew_adj %>%
      divide_cols(Xscale) %>%
      Matrix::Matrix(sparse = TRUE)
  }


  if (fit$backend == "e1071") {
    fit_svm <- fit$fit$best.model

    predict_obj <- predict(
      fit_svm,
      newdata = Xnew_adj,
      probability = TRUE
    )

    predict_prob <- attr(predict_obj, "probabilities") %>%
      as.matrix() %>%
      magrittr::set_rownames(
        rownames(Xnew_adj)
      ) %>%
      .[, sort(colnames(.))]

  } else if (fit$backend == "liquidSVM") {
    fit_svm <- fit$fit

    predict_prob <- predict(
      fit_svm,
      newdata = as.matrix(Xnew_adj)
      # probability = TRUE
    ) %>%
      data.matrix() %>%
      divide_rows(
        rowSums(.)
      ) %>%
      magrittr::set_colnames(
        colnames(.) %>%
          strsplit("vs") %>%
          sapply(head, 1)
      ) %>%
      magrittr::set_rownames(
        rownames(Xnew_adj)
      ) %>%
      .[, sort(colnames(.))]

  }


  pred_class <- apply(
    predict_prob,
    1,
    function(x) names(x)[which.max(x)]
  )

  list("predicted" = pred_class,
       "probs_predicted" = predict_prob,
       "observed" = Ynew)
}

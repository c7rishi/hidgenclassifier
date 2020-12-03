get_rand_foldid <- function(response, nfold = 10) {
  data.table::data.table(
    resp = response
  )[,
    foldid := sample(rep(1:nfold, length.out = .N)),
    by = resp
  ]$foldid
}


#' Hidden genome sparse multinomial logistic classifier (smlc)
#' @param X data design matrix with observations across rows and predictors across
#' columns. For a typical hidden genome classifier each row represents a tumor and
#' the columns represent (possibly normalized by some functions of
#' the total mutation burden in tumors) binary 1-0 presence/absence indicators
#' of raw variants, counts of mutations at specific genes and counts of mutations
#' corresponding to specific mutation signatures etc.
#' @param Y character vector or factor denoting the cancer type of tumors whose
#' mutation profiles are listed across the rows of \code{X}.
#' @param grouped logical. Use group-lasso penalty instead of the ordinary lasso
#' penalty? Defaults to TRUE.
#' @param ... additional arguments passed to \code{cv.glmnet}.
#'
#' @note
#' The function is a light wrapper around cv.glmnet with
#' \code{family = "multinomial"}, and \code{type.multinomial = "grouped"} if
#' \code{grouped} = TRUE. \code{cv.glmnet} tunes the sparsity hyper-parameter using
#' cross-validation. \code{fit_smlc} by default uses a 10-fold cross-validation
#' similar to the default of \code{cv.glmnet} (can be changed by supplying
#' \code{nfolds} in \code{...}); however with a stratified  random partition
#' (based on the categories of \code{Y}), instead of the default simple random
#' partition used in \code{cv.glmnet}. Override this by supplying \code{foldid} to
#' \code{cv.glmnet} in the \code{...}. In addition, \code{fit_smlc}
#' sets \code{maxit = 1e6}, \code{trace.it = TRUE} in \code{...} by default
#' (instead of the default
#' \code{maxit = 1e5} set in glmnet).
#'
#' @return
#' Returns a list containing the cv.glmnet fitted object,
#' the original X and Y and the estimated
#' intercept vector alpha and regression coefficients matrix beta.
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
#' # -- no meta-feature information used
#' fit0 <- fit_mlogit(
#'   X = var_design[idx_train, ],
#'   Y = canc_resp[idx_train]
#' )
#'
#' pred0 <- predict_mlogit(
#'   fit = fit0,
#'   Xnew = var_design[idx_test, ]
#' )
#'
#'
#' @export
fit_smlc <- function(X,
                     Y,
                     grouped = TRUE,
                     alpha = 1,
                     ...) {
  type.multinomial <- ifelse(grouped,
                             "grouped",
                             "ungrouped")
  dots <- list(...)
  dots$alpha <- NULL
  dots$type.multinomial <- NULL

  if (is.null(dots$nfolds)) {
    dots$nfolds <- 10
  }

  if (is.null(dots$foldid)) {
    dots$foldid <- get_rand_foldid(Y, dots$nfolds)
  }

  if (is.null(dots$maxit)) {
    dots$maxit <- 1e6
  }

  if (is.null(dots$keep)) {
    dots$keep <- TRUE
  }

  if (is.null(dots$trace.it)) {
    dots$trace.it <- TRUE
    if (!is.null(dots$parallel)) {
      if (dots$parallel) {
        dots$trace.it <- FALSE
      }
    }

  }

  logis <- do.call(
    glmnet::cv.glmnet,
    c(
      list(
        x = Matrix::Matrix(X, sparse = TRUE),
        y = Y,
        family = "multinomial",
        alpha = alpha,
        type.multinomial = type.multinomial
      ),
      dots
    )
  )


  tmp <- coef(logis)

  alpha_vec <- vapply(tmp, function(x) x["(Intercept)", ], 0)
  beta_mat <- do.call(cbind,
                      lapply(tmp, function(x) x[colnames(X), ]))


  list(alpha = alpha_vec,
       beta = beta_mat,
       X = X,
       Y = Y,
       fit = logis,
       method = "mlogit",
       glmnet_keep = dots$keep)
}

#' @rdname fit_smlc
#' @export
fit_mlogit <- fit_smlc

#' adjust Xnew by discarding columns not in Xold_colnames
#' and by adding 0-valued columns that are in Xold_colnames
#' but not in Xnew
adjust_Xnew <- function(Xnew, Xold_colnames) {
  all_preds <- union(colnames(Xnew), Xold_colnames)

  Xnew_adj <- Matrix::Matrix(
    0, nrow = nrow(Xnew),
    ncol = length(all_preds),
    dimnames = list(rownames(Xnew),
                    all_preds),
    sparse = TRUE)
  if (!is.null(rownames(Xnew))) {
    Xnew_adj[rownames(Xnew), colnames(Xnew)] <- Xnew
  } else {
    Xnew_adj[, colnames(Xnew)] <- Xnew
  }

  Xnew_adj <- Matrix::Matrix(Xnew_adj, sparse = TRUE)

  Xnew_adj
}


#' Prediction based on the hidden genome sparse multinomial
#' logistic classifier
#'
#' @param  Xnew test data design (or meta-design) matrix (observations
#' across rows and variables predictors/features across columns)
#' for which predictions are to be made from a fitted model. For a typical hidden
#' genome classifier this will be a matrix whose rows correspond to the test set
#' tumors, and columns correspond to (normalized by some functions of
#' the total mutation burdens in tumors) binary 1-0 presence/absence of
#' raw variants, counts of mutations at specific genes and counts of mutations
#' corresponding to specific mutation signatures etc.
#' @param fit fitted hidden genome classifier, an output of fit_smlc.
#' @param Ynew the actual cancer categories for the test samples.
#'
#' @note  Predictors in \code{Xnew} that are not present in the
#' training set design matrix (stored in \code{fit}) are dropped, and predictors
#' not included in \code{Xnew} but present in training set design matrix are
#' all assumed to have zero values. This is convenient for a typical
#' hidden genome classifier where most predictors are (some normalized versions
#' of) counts (e.g. for gene and mutation signatures) or
#' binary presence/absence indicators (e.g., for raw variants) so that a zero
#' predictor value essentially indicates some form of "absence".
#' However, care must be taken for predictors whose 0 values
#' do not indicate absence.
#' @seealso fit_mlogit
#'
#' @export
predict_smlc <- function(fit,
                         Xnew,
                         Ynew = NULL,
                         return_lin_pred = FALSE,
                         ...)  {
  fit_glmnet <- fit$fit
  beta <- fit$beta
  alpha <- fit$alpha

  Xnew_adj <- adjust_Xnew(Xnew, rownames(beta))
  n_new <- nrow(Xnew_adj)

  Xbeta_new <- (Xnew_adj[, rownames(beta), drop = FALSE] %*%
                  Matrix::Matrix(beta, sparse = TRUE) +
                  tcrossprod(rep(1, n_new), alpha))


  predict_prob <- t(apply(Xbeta_new, 1, softmax))


  out <- apply(predict_prob, 1, function(x) names(x)[which.max(x)])


  res <- list("predicted" = out,
              "probs_predicted" = predict_prob,
              "observed" = Ynew)
  if (return_lin_pred) {
    res$lin_pred <- Xbeta_new
  }

  res
}

#' @rdname predict_smlc
#' @export
predict_mlogit <- predict_smlc

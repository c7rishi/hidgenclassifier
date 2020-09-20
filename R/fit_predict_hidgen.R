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
#' @note   The function is a light wrapper around cv.glmnet with
#' \code{family = "multinomial"}, and \code{type.multinomial = "grouped"} if
#' \code{grouped} = TRUE.
#'
#'  @return
#' Returns a list containing the cv.glmnet fitted object,
#' the original X and Y and the estimated
#' intercept vector alpha and regression coefficients matrix beta.
#'
#' @export
fit_smlc <- function(X,  Y,
                     grouped = TRUE,
                     alpha = 1,
                     ...) {
  X <- Matrix::Matrix(X, sparse = TRUE)
  type.multinomial <- ifelse(grouped,
                             "grouped",
                             "ungrouped")
  dots <- list(...)
  dots$alpha <- NULL
  dots$type.multinomial <- NULL


  logis <- do.call(
    glmnet::cv.glmnet,
    c(
      list(x = X, y = Y,
           family = "multinomial",
           alpha = alpha,
           type.multinomial = type.multinomial),
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
       fit = logis)
}



#' Light wrapper around randomForest to use in hidden
#' genome classification
#' @inheritParams fit_smlc
#' @export
fit_rfc <- function(X, Y, ...) {
  randomForest::randomForest(
    x = as.matrix(X),
    y = as.factor(Y),
    ...
  )
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
#'
#' @export
predict_smlc <- function(Xnew,
                         fit,
                         Ynew = NULL,
                         return_lin_pred = FALSE,
                         ...)  {
  # browser()
  fit_glmnet <- fit$fit
  beta <- fit$beta
  alpha <- fit$alpha

  all_preds <- union(colnames(Xnew), rownames(beta))

  Xnew_adj <- Matrix::Matrix(0, nrow = nrow(Xnew),
                             ncol = length(all_preds),
                             dimnames = list(rownames(Xnew),
                                             all_preds),
                             sparse = TRUE)
  if (!is.null(rownames(Xnew))) {
    Xnew_adj[rownames(Xnew), colnames(Xnew)] <- Xnew
  } else {
    Xnew_adj[, colnames(Xnew)] <- Xnew
  }

  n_new <- nrow(Xnew_adj)
  Xnew_adj <- Matrix::Matrix(Xnew_adj, sparse = TRUE)

  # if (impute) {
  #   alpha_impute <- fit$alpha_impute
  #   beta0_impute <- fit$beta0_impute
  #   omega_impute <- fit$omega_impute # the SBS cols, to be scaled up
  #   allpred_impute <- rbind(beta0_impute, omega_impute)
  #
  #
  #   XU_non_target_impute <- exp(
  #     Xnew_adj[, rownames(allpred_impute),
  #              drop = FALSE] %*%
  #       allpred_impute +
  #       tcrossprod(rep(1, n_new), alpha_impute)
  #   ) + 1
  #
  #   # browser()
  #
  #   # scale up by adding the imputations
  #   Xnew_adj[, rownames(omega_impute)] <-
  #     Xnew_adj[, rownames(omega_impute)] +
  #     XU_non_target_impute[, rownames(omega_impute)]
  # }




  # X_common <- intersect(colnames(Xnew), rownames(beta))
  # X_absent <- setdiff(rownames(beta), colnames(Xnew))
  # if (length(X_absent) > 0) {
  #   zeromat <- matrix(0, nrow(Xnew),
  #                     length(X_absent)) %>%
  #     Matrix::Matrix(sparse = TRUE)
  #
  #   rownames(zeromat) <- rownames(Xnew)
  #   colnames(zeromat) <- X_absent
  #
  #   Xnew_adj <- cbind(Xnew[, X_common], zeromat)
  # } else {
  #   Xnew_adj <- Xnew[, X_common]
  # }
  #
  Xbeta_new <- (Xnew_adj[, rownames(beta), drop = FALSE] %*%
                  Matrix::Matrix(beta, sparse = TRUE) +
                  tcrossprod(rep(1, n_new), alpha))


  predict_prob <- t(apply(Xbeta_new, 1, softmax))

  # predict_prob <- predict(fit_glmnet,
  #                         newx = Xnew_adj[, rownames(beta)],
  #                         type = "response")[, , 1] %>%
  # as.matrix()
  # out1 <- apply(predict_prob1, 1, function(x) names(x)[which.max(x)])
  out <- apply(predict_prob, 1, function(x) names(x)[which.max(x)])


  res <- list("predicted" = out,
              "probs_predicted" = predict_prob,
              "observed" = Ynew)
  if (return_lin_pred) {
    res$lin_pred <- Xbeta_new
  }

  res
}



#
# fit_smlc_impute <- function(X, XU_target, XU_nontarget,
#                             Y,
#                             gaussian_alpha = 1,
#                             ...) {
#   `%>%` <- magrittr::`%>%`
#   # impute regression:
#   Xmat <- Matrix::Matrix(cbind(X, XU_target), sparse = TRUE)
#   Ymat <- log(1+as.matrix(XU_nontarget))
#
#   cat("\n Training imputation..")
#   impute_fit <- glmnet::cv.glmnet(
#     x = Xmat,
#     y = Ymat,
#     family = "mgaussian",
#     alpha = gaussian_alpha
#   )
#
#   param <- coef(impute_fit) %>%
#     do.call(cbind, .)
#
#   alpha_impute <- param["(Intercept)", ]
#   beta0_impute <- param[colnames(X), ]
#   omega_impute <- param[colnames(XU_target), ]
#   names(alpha_impute) <-
#     colnames(beta0_impute) <-
#     colnames(omega_impute) <-
#     colnames(Ymat)
#
#   cat("\n Training Classifier..")
#
#
#   res_smlc <- fit_smlc(
#     X = cbind(X,
#               (XU_target +
#                  XU_nontarget[, colnames(XU_target)])),
#     Y = Y, ...)
#
#
#   c(res_smlc,
#     list(
#       alpha_impute = alpha_impute,
#       beta0_impute = beta0_impute,
#       omega_impute = omega_impute,
#       X_impute = X,
#       XU_target_impute = XU_target,
#       XU_nontarget_impute = XU_nontarget
#     )
#   )
#
# }

# fit_lm_SBS <- function(Ymat, Zmat, gmat = NULL,
#                        ...) {
#   `%>%` <- magrittr::`%>%`
#   # impute regression:
#   Xmat <- log(1 + cbind(Zmat, gmat)) %>%
#     Matrix::Matrix(sparse = TRUE)
#
#   # browser()
#
#   logYmat <- log(1+as.matrix(Ymat))
#
#   # cat("\n Training non-target SBS using penalized lm..")
#   impute_fit <- glmnet::cv.glmnet(
#     x = Xmat,
#     y = logYmat,
#     family = "mgaussian",
#     ...
#   )
#
#   list(
#     fit = impute_fit,
#     Xmat = Xmat,
#     logYmat = logYmat
#   )
#
# }

# predict_lm_SBS <- function(lm_fit, Ymat_new, Zmat_new, gmat_new)



# .knn_predict_1obs <- function(xtest_single,
#                               Xmat_train,
#                               Ymat_train,
#                               k = sqrt(nrow(Xmat_train)),
#                               ...) {
#   k <- floor(k)
#   xtest_adj <- numeric(ncol(Xmat_train))
#   names(xtest_adj) <- colnames(Xmat_train)
#   xtest_adj[names(xtest_single)] <- xtest_single
#   dists_all <- apply(
#     Xmat_train,
#     1,
#     function(x) {
#       sqrt(sum((x - xtest_adj)^2))
#     }
#   )
#
#   rank_dist <- rank(dists_all)
#   knns <- which(rank_dist <= k)
#   colMeans(as.matrix(Ymat_train)[unname(knns), , drop = FALSE])
# }
#
#
#
# knn_prdict_SBS <- function(Xmat_train,
#                            Ymat_train,
#                            Xmat_test,
#                            Ymat_test = NULL,
#                            k = sqrt(nrow(Xmat_train)), ...) {
#   Ymat_predict <- t(
#     apply(
#       Xmat_test,
#       1,
#       function(x) {
#         .knn_predict_1obs(x, Xmat_train, Ymat_train, k = k)
#       }
#     )
#   )
#
#   if (!is.null(Ymat_test)) {
#     error <- Ymat_predict - Ymat_test
#
#     attr(Ymat_predict, "error") <- error
#   }
#
#   Ymat_predict
# }
#
#
# lm_predict_SBS <- function(fit, Ymat_new = NULL,
#                            Zmat_new, gmat_new = NULL, ...) {
#   `%>%` <- magrittr::`%>%`
#   if (!is.null(Ymat_new)) {
#     logYmat_new <- log(1+as.matrix(Ymat_new))
#   }
#   Xmat_new <- cbind(log(1+Zmat_new), gmat_new) %>%
#     Matrix::Matrix(sparse = TRUE)
#
#   predict_Y <- predict(fit$fit, newx = Xmat_new) %>%
#     .[, , 1] %>%
#     exp(.) - 1
#
#   if (!is.null(Ymat_new)) {
#     attr(predict_Y, "error") = predict_Y - Ymat_new
#     attr(predict_Y, "observed") = Ymat_new
#   }
#
#   predict_Y
#
# }


#' prediction based on hidden genome random forest classifier
#' @export
predict_rfc <- function(Xnew,
                        fit,
                        Ynew = NULL, ...)  {
  fit_rf <- fit
  X_common <- intersect(colnames(Xnew), rownames(fit_rf$importance))
  X_absent <- setdiff(rownames(fit_rf$importance), colnames(Xnew))
  if (length(X_absent) > 0) {
    zeromat <- matrix(0, nrow(Xnew),
                      length(X_absent))

    rownames(zeromat) <- rownames(Xnew)
    colnames(zeromat) <- X_absent

    Xnew_adj <- cbind(Xnew[, X_common], zeromat)
  } else {
    Xnew_adj <- Xnew[, X_common]
  }


  predict_prob <- as.matrix(
    predict(
      fit_rf,
      newdata = Xnew_adj[, rownames(fit_rf$importance)],
      type = "prob"
    )
  )[names(Ynew), ]

  out <- apply(predict_prob, 1, function(x) names(x)[which.max(x)])


  list("predicted" = out,
       "probs_predicted" = predict_prob,
       "observed" = Ynew)
}


#
# predict_smlc_component <- function(Xnew,
#                                    fit,
#                                    Ynew = NULL,
#                                    impute = FALSE,
#                                    impute_pv = NULL,
#                                    impute_pvU = NULL,
#                                    ...)  {
#   # browser()
#   `%>%` <- magrittr::`%>%`
#   fit_glmnet <- fit$fit
#   beta <- fit$beta
#   alpha <- fit$alpha
#
#   all_preds <- union(colnames(Xnew), rownames(beta))
#
#   Xnew_adj <- Matrix::Matrix(0, nrow = nrow(Xnew),
#                              ncol = length(all_preds),
#                              dimnames = list(rownames(Xnew),
#                                              all_preds),
#                              sparse = TRUE)
#
#   Xnew_adj[rownames(Xnew), colnames(Xnew)] <- Xnew
#   n_new <- nrow(Xnew_adj)
#
#   if (impute) {
#     alpha_impute <- fit$alpha_impute
#     beta0_impute <- fit$beta0_impute
#     omega_impute <- fit$omega_impute # the SBS cols, to be scaled up
#     allpred_impute <- rbind(beta0_impute, omega_impute)
#
#
#     XU_non_target_impute <- exp(
#       Xnew_adj[, rownames(allpred_impute),
#                drop = FALSE] %*%
#         allpred_impute +
#         tcrossprod(rep(1, n_new), alpha_impute)
#     ) + 1
#
#     browser()
#
#     # scale up by adding the imputations
#     Xnew_adj[, rownames(omega_impute)] <-
#       Xnew_adj[, rownames(omega_impute)] +
#       XU_non_target_impute[, rownames(omega_impute)]
#   }
#
#
#
#   # browser()
#
#
#   Xbeta_comp <- lapply(
#     1:nrow(Xnew_adj),
#     function(j) {
#       rbind(
#         intercept = alpha,
#         beta * Xnew_adj[j, rownames(beta)]
#       )
#       # same as
#       # apply(beta, 2, function(x) x * Xnew_adj[j, rownames(beta)])
#     }
#   ) %>%
#     setNames(rownames(Xnew_adj))
#
#   # Xbeta_new <- (Xnew_adj[, rownames(beta), drop = FALSE] %*% beta +
#   #                 tcrossprod(rep(1, n_new), alpha))
#   #
#   #
#   # predict_prob <- t(apply(Xbeta_new, 1, softmax))
#
#   # predict_prob <- predict(fit_glmnet,
#   #                         newx = Xnew_adj[, rownames(beta)],
#   #                         type = "response")[, , 1] %>%
#   # as.matrix()
#   # out1 <- apply(predict_prob1, 1, function(x) names(x)[which.max(x)])
#   # out <- apply(predict_prob, 1, function(x) names(x)[which.max(x)])
#   #
#   #
#   # list("predicted" = out,
#   #      "probs_predicted" = predict_prob,
#   #      "observed" = Ynew)
#   Xbeta_comp
# }

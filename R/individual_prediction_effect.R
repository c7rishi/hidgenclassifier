# individual component terms in smlc_prediction,
# used for determining prediction effects of each
# variable in *each* test set tumor
predict_wo_feature <- function(fit,
                               Xnew,
                               ...)  {
  if (fit$method %in% c("mlogit")) {
    fit_glmnet <- fit$fit
    beta <- fit$beta
    alpha <- fit$alpha

    Xnew_adj <- adjust_Xnew(Xnew, rownames(beta))
    n_new <- nrow(Xnew_adj)

    Xbeta_new <- (Xnew_adj[, rownames(beta), drop = FALSE] %*%
                    Matrix::Matrix(beta, sparse = TRUE) +
                    tcrossprod(rep(1, n_new), alpha))

    Xbeta_comp_list <- lapply(
      1:nrow(Xnew_adj),
      function(j) {
        rbind(
          intercept = alpha,
          beta * Xnew_adj[j, rownames(beta)]
        ) %>%
          .[!grepl("intercept", rownames(.), ignore.case = TRUE), ]
        # same as
        # apply(beta, 2, function(x) x * Xnew_adj[j, rownames(beta)])
      }
    ) %>%
      setNames(rownames(Xnew_adj))

    col_sums_list <- Xbeta_comp_list %>%
      lapply(function(x) apply(x, 2, sum))
    full_pred_list <- lapply(col_sums_list, softmax)

    pred_wo_variable <- mapply(
      function(col_sums, Xbeta_comp) {
        (tcrossprod(rep(1, nrow(Xbeta_comp)), col_sums) - Xbeta_comp) %>%
          apply(1, softmax) %>%
          t()
      },
      col_sums = col_sums_list,
      Xbeta_comp = Xbeta_comp_list,
      SIMPLIFY = FALSE
    )
  }
  list(pred_wo_variable = pred_wo_variable,
       full_pred = full_pred_list)
}


# kld <- function(p, q) {
#   p_non0 <- which(p > 0)
#   pmax(sum(p * log(pmax(p, 1e-10)/pmax(q, 1e-10))), 0)
# }
#
# jsdist <- function(p, q) {
#   m <- (p + q)/2
#   js <- pmax((kld(p, m) + kld(q, m))/2, 0)
#   sqrt(js/log(2))
# }


#' Probability distance based effects of predictors in a fitted model
#' in each *individual prediction*
#'
#' @inheritParams fit_smlc
#' @param dist Probability distance measure to use. Defaults to "jsdist"
#' which is the Jensen-Shannon distance
#'
#' @export
indiv_predict_effect <- function(fit, Xnew, dist = "jsdist", ...) {
  # evaluate classification probabilities for Xnew when
  # each predictor is set to zero
  if (!fit$method %in% c("mlogit")) {
    stop("Only implemented for mlogit classifiers yet")
  }
  pred_prob_list <- predict_wo_feature(fit = fit, Xnew = Xnew)
  outmat <- mapply(
    function(pred_wo_variable, full_pred) {
      calc_dist_refvec_targetmat(P = full_pred, Qmat = pred_wo_variable) %>%
        setNames(rownames(pred_wo_variable))
    },
    pred_wo_variable = pred_prob_list$pred_wo_variable,
    full_pred = pred_prob_list$full_pred,
    SIMPLIFY = FALSE
  ) %>%
    do.call(rbind, .) %>%
    Matrix::Matrix(sparse = TRUE)

  outmat
}


# row_contrib_to_colsums_softmax <- function(mat) {
#   col_sums <- apply(mat, 2, sum)
#   full_pred <- softmax(col_sums)
#   # p = pred_wo_row["cyto_chr_18_cyto_p11.31", ]
#   # jsdist(p, full_pred)
#   # browser()
#   pred_wo_row <- (tcrossprod(rep(1, nrow(mat)), col_sums) - mat) %>%
#     apply(1, softmax) %>%
#     t()
#   jsdist_vals <- apply(pred_wo_row, 1, jsdist, q = full_pred)
#
#   cbind(pred_wo_row, "JS_dist" = jsdist_vals) %>%
#     tibble::as_tibble(rownames = "Predictors") %>%
#     dplyr::arrange(dplyr::desc(JS_dist))
#   # dplyr::filter(!grepl("intercept", Predictors, ignore.case = TRUE))
# }

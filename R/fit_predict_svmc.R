#' Hidden genome SVM classifier (svmc)
#'
#' @details Light wrapper around e1071::svm to use in hidden
#' genome classification
#' @param ... additional arguments passed to e1071:tune.svm.
#' @inheritParams fit_smlc
#' @export
fit_svmc <- function(X, Y, ...) {
  dots <- list(...)
  if (is.null(dots$gamma)) {
    gamma <- 10^(-4:3)
  }
  if (is.null(dots$cost)) {
    dots$cost <- 10^(-4:3)
  }

  fit <- list(
    x = X,
    y = as.factor(Y),
    probability = TRUE,
    kernel = "radial"
  ) %>%
    c(dots) %>%
    do.call(e1071::tune.svm, .)

  out <- list(
    X = X,
    Y = Y,
    fit = fit
  )
}




#' prediction based on hidden genome random forest classifier
#' @export
predict_svmc <- function(Xnew,
                         fit,
                         Ynew = NULL, ...)  {

  fit_svm <- fit$fit$best.model

  Xold_names <- colnames(fit$X)
  Xnew_adj <- Xnew %>%
    fill_sparsemat_zero(
      rownames = rownames(.),
      colnames = Xold_names
    )

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


  pred_class <- apply(
    predict_prob,
    1,
    function(x) names(x)[which.max(x)]
  )

  list("predicted" = pred_class,
       "probs_predicted" = predict_prob,
       "observed" = Ynew)
}

#' Hidden genome random forest classifier (rfc)
#'
#' @details Light wrapper around randomForest or ranger to use in hidden
#' genome classification
#' @param backend Which backend to use? Available options are
#' "ranger" and "randomForest" corresponding to the respective R packages.
#' NOTE: randomForest does not support sparseMatrix, and the predictor matrix
#' is coerced into an ordinary matrix. This means using randomForest will likely
#' be more memory intensive and hence slower than ranger.
#' @param ... additional arguments passed to ranger or randomForest (depending
#' on backend).
#' @inheritParams fit_smlc
#' @export
fit_rfc <- function(X, Y, backend = "ranger", ...) {

  if (!backend %in% c("randomForest", "ranger")) {
    stop('backend must be one of "randomForest" or "ranger"')
  }

  if (backend == "ranger") {
    fit <- ranger::ranger(
      x = X,
      y = as.factor(Y),
      classification = TRUE,
      probability = TRUE,
      ...
    )
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
    backend = backend
  )

  out
}




#' prediction based on hidden genome random forest classifier
#' @export
predict_rfc <- function(fit,
                        Xnew,
                        Ynew = NULL, ...)  {

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

  list("predicted" = pred_class,
       "probs_predicted" = predict_prob,
       "observed" = Ynew)
}

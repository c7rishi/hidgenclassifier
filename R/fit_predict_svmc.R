#' Hidden genome SVM classifier (svmc)
#'
#' @details Light wrapper around e1071::svm or liquidSVM::mcSVM to use in hidden
#' genome classification
#' @param ... additional arguments passed to e1071:tune.svm, or
#' liquidSVM::svm.
#' @param backend the backend to use. Either "e1071" or "liquidSVM". Defaults to
#' "liquidSVM"
#' @inheritParams fit_smlc
#' @export
fit_svmc <- function(X,
                     Y,
                     backend = "liquidSVM",
                     scale = TRUE,
                     ...) {
  dots <- list(...)
  if (is.null(dots$gamma)) {
    gamma <- 10^(-4:3)
  }
  if (is.null(dots$cost)) {
    dots$cost <- 10^(-4:3)
  }


  if (scale) {
    X_scale <- apply(X, 2, max) %>%
      pmax(1) %>%
      # ifelse(. < 1e-10, 1, .) %>%
      setNames(colnames(X))
    X <- X %>%
      scale(center = FALSE, scale = X_scale) %>%
      Matrix::Matrix(sparse = TRUE)

    # X[is.nan(X)] <- 0
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
      dots$min_gamma <- 1e-4
    }

    if (is.null(dots$min_lambda)) {
      dots$min_lambda <- 1e-7
    }

    # if (is.null(dots$scale)) {
    #   dots$scale <- FALSE
    # }


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
    backend = backend
  )

  out
}




#' prediction based on hidden genome random forest classifier
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

    Xscale[colnames(Xold_names)] <- attr(
      fit$X,
      "scaled:scale"
    )[colnames(fit$X)]

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

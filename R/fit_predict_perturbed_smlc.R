# Perturbed fit of penalized multilogit
# based on an existing fitted penalized multilogit
fit_perturbed_smlc <- function(fit,
                               exact = FALSE,
                               lambda_type = "lambda.1se",
                               ...) {
  dots <- list(...)
  dots$alpha <- NULL
  dots$type.multinomial <- NULL

  glment_alpha <- fit$glmnet_alpha
  if (is.null(glment_alpha)) glmnet_alpha <- 1

  glment_type.multinomial <- fit$type.multinomial
  if (is.null(glment_type.multinomial)) glmnet_type.multinomial <- "grouped"

  stopifnot(
    lambda_type %in% c("lambda.min", "lambda.1se"),
    length(lambda_type) == 1,
    is.logical(exact),
    length(exact) == 1
  )
  lambda <- lambda_orig <- fit$fit[[lambda_type]]
  alpha <- fit$alpha

  if (is.null(dots$maxit)) {
    dots$maxit <- 1e6
  }

  if (is.null(dots$trace.it)) {
    dots$trace.it <- TRUE
    if (!is.null(dots$parallel)) {
      if (dots$parallel) {
        dots$trace.it <- FALSE
      }
    }
  }


  wt_lambda <- 1 # rexp(1)
  wt_data <- rexp(nrow(fit$X))
  wt_prior <- rexp(ncol(fit$X))
  lambda_final <- lambda_orig * wt_lambda

  dots$lambda <- NULL
  dots$weights <- wt_data
  dots$penalty.factor <- wt_prior

  perturb_weights <- list(
    data = wt_data,
    prior = wt_prior,
    lambda = wt_lambda
  )

  logis <- do.call(
    glmnet::glmnet,
    c(
      list(
        x = Matrix::Matrix(fit$X, sparse = TRUE),
        y = fit$Y,
        family = "multinomial",
        alpha = glmnet_alpha,
        type.multinomial = glmnet_type.multinomial
      ),
      dots
    )
  )


  tmp <- coef(
    logis,
    s = lambda_final,
    exact = exact,
    x = Matrix::Matrix(fit$X, sparse = TRUE),
    y = fit$Y,
    weights = wt_data,
    penalty.factor = wt_prior
  )

  icept_idx <- 1
  alpha_vec <- sapply(tmp, function(xx) xx[icept_idx, ])

  beta_mat <- do.call(
    cbind,
    lapply(tmp, function(xx) xx[colnames(fit$X), ])
  )


  res <- list(
    alpha = alpha_vec,
    beta = beta_mat,
    X = fit$X,
    Y = fit$Y,
    lambda = lambda_orig,
    fit = logis,
    method = "mlogit",
    glmnet_keep = dots$keep,
    perturb_weights = perturb_weights,
    glmnet_alpha = glmnet_alpha,
    glmnet_type.multinomial = glmnet_type.multinomial
  )

  X_sc <- scale(fit$X, center = FALSE, scale = TRUE)
  X_col_sd <- attr(X_sc, "scaled:scale")

  lpd_both <- list(
    orig = fit,
    curr = res
  ) %>%
    sapply(
      function(this_fit) {
        calc_lpd_smlc(
          fit = this_fit,
          lambda = lambda_final,
          X_col_sd = X_col_sd,
          glmnet_alpha = glmnet_alpha,
          glmnet_type.multinomial = glmnet_type.multinomial
        )
      }
    )
  res$lpd <- lpd_both
  res$log_importance <- unname(lpd_both[1] - lpd_both[2])

  res
}


calc_lpd_smlc <- function(fit,
                          lambda,
                          X_col_sd = NULL,
                          glmnet_alpha,
                          glmnet_type.multinomial,
                          ...) {
  glment_predict_prob <- predict(
    fit$fit,
    newx = Matrix::Matrix(fit$X, sparse = TRUE),
    type = "response"
  )[, , 1]

  N <- nrow(fit$X)
  loglik <- sapply(
    1:N,
    function(ii) {
      glment_predict_prob[ii, unname(fit$Y[ii])]
    }
  ) %>%
    log(.) %>%
    sum(.)

  if (is.null(X_col_sd)) {
    X_sc <- scale(fit$X, center = FALSE, scale = TRUE)
    X_col_sd <- attr(X_sc, "scaled:scale")
  }

  # need to rescale rows of beta by X_col_sd
  beta <- fit$beta * X_col_sd
  # equivalent to fit$beta %>% divide_rows(1/X_col_sd)
  norm_beta_F <- sqrt(sum(beta^2))

  norm_beta_group <- ifelse(
    glmnet_type.multinomial == "grouped",
    apply(beta, 1, function(xx) sqrt(sum(xx^2))) %>% sum(),
    sum(abs(beta))
  )

  # need to multiply with the sample size N as glmnet
  # uses a normalized version of the objective funciton
  lprior <- - lambda * N * (
    (1 - glmnet_alpha) * norm_beta_F / 2 +
      glmnet_alpha * norm_beta_group
  )

  loglik + lprior
}

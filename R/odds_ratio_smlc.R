#' Calculate odds ratios from a multinomial logistic
#' hidden genome model
#' @inheritParams predict_mlogit
#' @param type either "one-vs-rest" or "one-vs-one".
#' @param baseline_category The response category with respect to which
#' odds ratios are calculated from the fitted model if \code{type = "one-vs-one"}.
#' If NULL (default), the first response cancer category obtained after sorting the
#' column labels (using \code{sort()}) is used.
#' @param log logical. Should the odds ratios be returned in log scale?
#' Defaults to TRUE.
#' @param predictor_subset Character vector listing the subset of predictors in
#' the fitted model \code{fit} for which odds ratios are to be computed.
#' @details
#' If  \code{type = "one-vs-rest"}, odds ratios are calculated for each predictor,
#' across all response categories, for one unit (one standard deviation unit,
#' if \code{scale = TRUE}) increase in each predictor at its mean, while keeping
#' all other predictors fixed at their respective means.
#'
#' If \code{type = "one-vs-one"}, odds ratios relative to a baseline category
#' is calculated. (Not implemented yet.)
#'
#' @return
#' Returns a sparse matrix (of class dgeMatrix) with odds ratios for predictors
#' (along the rows) across cancer sites (along the columns).
#'
#' @export
odds_ratio_mlogit <- function(
  fit,
  type = "one-vs-rest",
  scale = TRUE,
  log = TRUE,
  predictor_subset = NULL,
  baseline_category = NULL,
  ...
) {

  if (is.null(predictor_subset)) {
    predictor_subset <- colnames(fit$X)
  }

  Xmat <- fit$X
  Xmat_scale <- scale(Xmat)

  all_type <- c("one-vs-rest")

  if (!type %in% all_type) {
    msg <- paste("'type' must be one of",
                 paste(all_type, collapse = ", "))
    stop(msg)
  }

  if (type == "one-vs-rest") {
    mu <- attr(Xmat_scale, "scaled:center")
    sigma <- attr(Xmat_scale, "scaled:scale")
    d <- length(mu)
    d1 <- length(predictor_subset)

    x0 <- tcrossprod(rep(1, d), mu)
    x1 <- x0 + diag(sigma, nrow = d)
    dimnames(x0) <- dimnames(x1) <- list(names(mu), names(mu))

    x1 <- x1[predictor_subset, ]
    x0 <- x0[predictor_subset, ]

    x1beta <- predict_smlc(
      Xnew = x1,
      fit = fit,
      return_lin_pred = TRUE
    )$lin_pred

    x0beta <- predict_smlc(
      Xnew = x0[1, , drop = FALSE],
      fit = fit,
      return_lin_pred = TRUE
    )$lin_pred %>%
      as.vector() %>%
      tcrossprod(
        rep(1, d1),
        .
      )
    dimnames(x0beta) <- dimnames(x1beta)

    adj_const <- max(x1beta, x0beta)

    term3 <- .log_exp_shift_sum_rest_cols(x1beta, shift = adj_const)
    term4 <- .log_exp_shift_sum_rest_cols(
      x0beta[1, , drop = FALSE],
      shift = adj_const
    ) %>%
      as.numeric() %>%
      tcrossprod(
        rep(1, d1),
        .
      )
    dimnames(term4) <- dimnames(term4)

    out <- x1beta - x0beta - term3 + term4
  }
  if (!log) {
    out <- exp(out)
  }

  out
}

#' @rdname odds_ratio_mlogit
#' @export
odds_ratio_smlc <- odds_ratio_mlogit

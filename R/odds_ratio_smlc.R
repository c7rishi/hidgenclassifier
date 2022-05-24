#' Calculate odds ratios from a multinomial logistic
#' hidden genome model
#' @inheritParams predict_mlogit
#' @param type either "one-vs-rest" (default) or "one-vs-ave-baseline".
#' @param baseline_category Vector of response categories with respect to whose
#' *geometric average probability* are the generalized odds calculated.
#' Ignored if `type = "one-vs-rest"`
#' @param exclude_itself_from_baseline logical. If `type = "one-vs-ave-baseline"`
#' should the response category category site whose odds is being calculated
#' (e.g., category `B` in the numerator in the formula given in Details)
#' be excluded from the categories specified
#' in `baseline_category`? Defaults to TRUE. Ignored if `type = "one-vs-rest"` or
#' `length(baseline_category) == 1`.
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
#' If \code{type = "one-vs-ave-baseline"}, ratio of generalized odds ratios
#' relative to the *geometric average* of baseline category probabilities
#' are computed. For example if `baseline_category = c("A1", .., "Ak")`
#' then the generalized odds of response
#' cancer site `B`is defined as $Pr(Site = B)/(\prod_{h=1}^k Pr(Site = Ak))^{1/k}$.
#' odds ratios are calculated from the fitted model if \code{type = "one-vs-one"}.
#' If NULL (default), all the response cancer categories are used.
#'
#' @return
#' Returns a sparse matrix (of class dgeMatrix) with odds ratios for predictors
#' (along the rows) across cancer sites (along the columns).
#'
#' @md
#'
#' @export
odds_ratio_mlogit <- function(
  fit,
  type = c("one-vs-rest", "one-vs-ave-baseline"),
  scale_1sd = TRUE,
  log = TRUE,
  predictor_subset = NULL,
  baseline_category = NULL,
  exclude_itself_from_baseline = TRUE,
  ...
) {

  if (is.null(predictor_subset)) {
    predictor_subset <- colnames(fit$X)
  }

  type <- match.arg(type)


  all_type <- c("one-vs-rest", "one-vs-ave-baseline")

  if (!type %in% all_type) {
    msg <- paste("'type' must be one of",
                 paste(all_type, collapse = ", "))
    stop(msg)
  }

  if (type == "one-vs-rest") {
    Xmat <- fit$X
    Xmat_scale <- scale(Xmat)

    mu <- attr(Xmat_scale, "scaled:center")
    sigma <- attr(Xmat_scale, "scaled:scale")
    d <- length(mu)
    d1 <- length(predictor_subset)

    x0 <- tcrossprod(rep(1, d), mu)
    x1 <- (x0) %>%
      `diag<-`(diag(.) + sigma)
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
    dimnames(term4) <- dimnames(term3)

    out <- x1beta - x0beta - term3 + term4
  } else if (type == "one-vs-ave-baseline") {

    if (is.null(baseline_category)) {
      baseline_category <- colnames(fit$beta)
    }

    all_cat <- colnames(fit$beta)

    no_match_cat <- setdiff(baseline_category, all_cat)

    if (length(no_match_cat) > 0) {
      msg <- paste0("'", no_match_cat, "'") %>%
        paste("baseline categories", .) %>%
        paste("not found")
      stop(msg)
    }

    stopifnot(length(baseline_category) > 0)

    betamat <- fit$beta

    if (scale_1sd) {
      fit1 <- fit
      fit$X <- scale(fit$X, center = FALSE)
      X_sc <- attr(fit$X, "scaled:scale") %>%
        {ifelse(. > 0, ., 1)}
      betamat <- betamat %>%
        divide_rows(1/c(X_sc[rownames(.)]))
    }


    exclude_itself <- exclude_itself_from_baseline &
      length(baseline_category > 1)

    out <- betamat
    for (jj in all_cat) {
      this_baseline <- baseline_category %>%
        {if (exclude_itself) setdiff(., jj) else .}
      out[, jj] <- out[, jj] - rowMeans(out[, this_baseline])
    }

  }


  if (!log) {
    out <- exp(out)
  }

  attr(out, "type") <- type

  out
}

#' @rdname odds_ratio_mlogit
#' @export
odds_ratio_smlc <- odds_ratio_mlogit

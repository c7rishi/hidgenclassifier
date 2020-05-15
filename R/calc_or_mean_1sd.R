
#' @export
calc_OR_mean_1SD <- function(fit, ...) {
  `%>%` <- magrittr::`%>%`
  Xmat <- fit$X
  Xmat_scale <- scale(Xmat)

  mu <- attr(Xmat_scale, "scaled:center")
  sigma <- attr(Xmat_scale, "scaled:scale")
  d <- length(mu)


  x0 <- tcrossprod(rep(1, d), mu)
  x1 <- x0 + diag(sigma, nrow = d)
  dimnames(x0) <- dimnames(x1) <- list(names(mu), names(mu))

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
      rep(1, d),
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
      rep(1, d),
      .
    )
  dimnames(term4) <- dimnames(term4)

  out <- x1beta - x0beta - term3 + term4

  out
}

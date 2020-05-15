# supporting functions
# random folds via stratified sampling
get_rand_foldid <- function(response, nfold = 10) {
  data.table::data.table(
    resp = response
  )[,
    foldid := sample(rep(1:nfold, length.out = .N)),
    by = resp
    ]$foldid
}


softmax <- function(x) {
  x1 <- exp(x - max(x))
  s_x1 <- sum(x1)
  if (s_x1 > 0) {
    x1/s_x1
  } else {
    x1
  }
}




.log_sum_exp <- function(x) {
  max_x <- max(x)
  if (max_x > 0) {
    log(max_x) + log(sum(exp(x-max_x)))
  } else {
    log(sum(exp(x)))
  }
}


.log_exp_shift_sum_rest_cols <- function(x, shift) {
  exp_x <- exp(x - shift)
  nc <- ncol(exp_x)
  rs_x <- rowSums(exp_x)
  exp_out <- tcrossprod(rs_x, rep(1, nc)) - exp_x
  dimnames(exp_out) <- dimnames(exp_x)
  log(exp_out)
}



harmonic_mean <- function(...) {
  dots <- list(...)
  x <- unlist(dots)
  1/mean(1/x)
}



kld <- function(p, q) {
  p_non0 <- which(p > 0)
  pmax(sum(p * log(pmax(p, 1e-10)/pmax(q, 1e-10))), 0)
}

jsdist <- function(p, q) {
  m <- (p + q)/2
  js <- pmax((kld(p, m) + kld(q, m))/2, 0)
  sqrt(js/log(2))
}



# row_contrib_to_colsums_softmax <- function(mat) {
#   col_sums <- apply(mat, 2, sum)
#   full_pred <- softmax(col_sums)
#
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

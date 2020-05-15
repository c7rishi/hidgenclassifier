#' calculates normalized mutual information given P(variant | cancer)
#' and p(cancer)
#' @export
calc_minfo <- function(prob_v_given_c, prob_c, normalize = TRUE) {
  prob_c <- prob_c/sum(prob_c)
  mi <- Cpp_calc_minfo(prob_v_given_c,
                       prob_c[colnames(prob_v_given_c)],
                       as.integer(normalize))
  names(mi) <- rownames(prob_v_given_c)
  mi
}

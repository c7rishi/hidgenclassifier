#include <Rcpp.h>
using namespace Rcpp;


double xlogx (double x) {
  if (x > 0) {
    return (x * log(x));
  } else {
    return 0.0;
  }
}



// returns 0 if x <= 0
double zero_log (double x) {
  if (x > 0) {
    return log(x);
  } else {
    return 0.0;
  }
}



// [[Rcpp::export]]
NumericVector Cpp_calc_minfo(NumericMatrix prob_mat,
                             NumericVector wt_vec,
                             int normalized = 1) {
  unsigned int n = prob_mat.nrow(), K = prob_mat.ncol();
  double p_joint_1, p_joint_0,
  p_row_marginal_1, p_row_marginal_0, tmp;

  NumericVector log_wt_vec = log(wt_vec);

  double sqrt_H_wt = 1.0, sqrt_H_row_marginal = 1.0;


  if (normalized == 1) {
    double tmp_sum = 0.0;
    for (unsigned int i = 0; i < K; i++) {
      tmp_sum -= log_wt_vec[i] * wt_vec[i];
    }

    sqrt_H_wt = sqrt(tmp_sum);
  }


  // Rcpp::Rcout << "sqrt_H_wt: " << sqrt_H_wt << std::endl;

  NumericVector minfo_vec(n);

  for(unsigned int i = 0; i < n; i++) {
    if (i % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    p_row_marginal_0 = 0;
    p_row_marginal_1 = 0;

    tmp = 0;

    for(unsigned int j = 0; j < K; j++) {
      p_joint_1 = prob_mat(i, j) * wt_vec[j];
      p_joint_0 = (1-prob_mat(i, j)) * wt_vec[j];

      p_row_marginal_1 += p_joint_1;
      p_row_marginal_0 += p_joint_0;

      tmp += xlogx(p_joint_1) - p_joint_1 * log_wt_vec[j] +
        xlogx(p_joint_0) - p_joint_0 * log_wt_vec[j];
    }

    minfo_vec[i] = tmp - xlogx(p_row_marginal_1) - xlogx(p_row_marginal_0);

    if (normalized == 1) {
      sqrt_H_row_marginal = sqrt(
        - (xlogx(p_row_marginal_0) + xlogx(p_row_marginal_1))
      );

      minfo_vec[i] /= (sqrt_H_row_marginal * sqrt_H_wt);
    }
    // Rcpp::Rcout << "NMI: " << minfo_vec[i] << std::endl;
  }


  return minfo_vec;
}

#include <Rcpp.h>
using namespace Rcpp;


double xlogyplus(double x, double y) {
  if (y > 0) {
    return (x * log(y));
  } else {
    return 0.0;
  }
}


double kldiv_e_pair(NumericVector P, NumericVector Q) {
  const int n = P.length();
  double kl = 0.0;
  for (int i = 0; i < n; i++) {
    kl += xlogyplus(P(i), P(i)) - xlogyplus(P(i), Q(i));
  }
  return kl;
}

double shannon(NumericVector P) {
  const int n = P.length();
  double kl = 0.0;
  for (int i = 0; i < n; i++) {
    kl -= xlogyplus(P(i), P(i));
  }
  return kl;
}

double jsdiv_pair(NumericVector P, NumericVector Q) {
  NumericVector M = (P + Q) / 2.0;
  double jsdiv = 0.5 * (kldiv_e_pair(P, M) + kldiv_e_pair(Q, M)) / log(2.0);
  return jsdiv;
}


double jsdiv_pair_normalized(NumericVector P, NumericVector Q) {
  NumericVector M = (P + Q) / 2.0;
  double jsdiv = 0.5 * (kldiv_e_pair(P, M) + kldiv_e_pair(Q, M)) / sqrt(shannon(M));
  return jsdiv;
}


double hellinger_pair(NumericVector P, NumericVector Q) {
  double sum = 0.0;
  const int n = P.length();
  for (int i = 0; i < n; i++) {
    sum += sqrt(P(i) * Q(i));
  }
  return sqrt(1 - sum);
}


double tv_pair(NumericVector P, NumericVector Q) {
  double sum = 0.0;
  const int n = P.length();
  for (int i = 0; i < n; i++) {
    sum += fabs(P(i) - Q(i));
  }
  return 0.5 * sum;
}


double renyi_entropy(NumericVector P, double alpha = 2.0) {
  if (fabs(alpha - 1.0) < 1e-20) {
    return shannon(P);
  } else {
    const int n = P.length();
    double tmp = 0.0;
    for (int i = 0; i < n; i++) {
      tmp += pow(P(i), alpha);
    }
    return log(tmp)/(1-alpha);
  }
}

double jensen_renyi_pair(NumericVector P, NumericVector Q, double alpha = 2) {
  NumericVector M = (P + Q) / 2.0;
  double out =
    renyi_entropy(M, alpha) -
    0.5 * (
        renyi_entropy(P, alpha) +
          renyi_entropy(Q, alpha)
    );

  return out;
}


// [[Rcpp::export]]
double calc_dist_pair(NumericVector P,
                      NumericVector Q,
                      String dist_type = "jsdist",
                      double alpha = 2) {
  double out = 0.0;

  if (dist_type == "jsdist") out = sqrt(jsdiv_pair(P, Q));
  else if (dist_type == "jsdiv") out = jsdiv_pair(P, Q);
  else if (dist_type == "jsdiv_normalized") out = jsdiv_pair_normalized(P, Q);
  else if (dist_type == "hellinger") out = hellinger_pair(P, Q);
  else if (dist_type == "tv") out = tv_pair(P, Q);
  else if (dist_type == "jrdiv") out = jensen_renyi_pair(P, Q, alpha);

  return out;
}

// P and each row of Qmat is a prob vector
// [[Rcpp::export]]
NumericVector calc_dist_refvec_targetmat(NumericVector P,
                                         NumericMatrix Qmat,
                                         String dist_type = "jsdist",
                                         double alpha = 2) {
  int n = Qmat.nrow(), k = Qmat.ncol();
  NumericVector out(n);
  NumericVector tmp(k);
  for (int i = 0; i < n; i++) {
    // extract ith row from Qmat
    for (int j = 0; j < k; j++) {
      tmp(j) = Qmat(i, j);
    }
    out(i) = calc_dist_pair(P, tmp, dist_type, alpha);
  }
  return(out);
}

// // [[Rcpp::export]]
// NumericMatrix calc_dist_mat(Rcpp::List Plist,
//                             String dist_type = "jsdist",
//                             double alpha = 1.0) {
//   const int n = Plist.size();
//   NumericMatrix dist = arma::zeros(n, n);
//   for (int i = 0; i < n - 1; i++) {
//     for (int j = i+1; j < n; j++) {
//
//       dist(i, j) = calc_dist_pair(Plist[i], Plist[j], dist_type, alpha);
//       dist(j, i) = dist(i, j);
//     }
//   }
//   return dist;
// }

## usethis namespace: start
#' @useDynLib hidgenclassifier, .registration = TRUE
#' @import data.table
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom Matrix Matrix crossprod tcrossprod t rowSums colSums rowMeans colMeans
#' @import stats
#' @import methods
#' @import keras
#' @importFrom ggplot2 autoplot
## usethis namespace: end
NULL

`%>%` <- magrittr::`%>%`

#' Convenience functions for dividing rows and columns of a matrix
#' by a given vector
#'
#' @param mat A matrix. Can be a sparseMatrix from {Matrix}
#' @param vec
#'
#' @details
#' \code{divide_rows} (\code{divide_cols}) divide the i-th row
#' (i-th column) of \code{mat} by the i-th entry of \code{vec}.
#' Formally, \code{divide_rows(mat, vec)} is equivalent to
#' sweep(mat, 1, vec "/"), but is usually more efficient, and
#' plays nicer with {magrittr} pipe \code{%>%}.
#'
#' @export
divide_rows <- function(mat, vec) {
  mat * (1/vec)
}

#' @rdname divide_rows
#' @export
divide_cols<- function(mat, vec) {
  t(divide_rows(t(mat), vec))
}

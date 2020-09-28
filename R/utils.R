# expand a sparse matrix by filling zeros, given rows & colnames
fill_sparsemat_zero <- function(
  mat, rownames, colnames
) {
  # mat_sparse <- Matrix::Matrix(mat, sparse = TRUE)
  zeromat <- Matrix(
    data = 0,
    nrow = length(rownames),
    ncol = length(colnames),
    dimnames = list(rownames, colnames),
    sparse = TRUE
  )
  rownames_common <- intersect(rownames, rownames(mat))
  colnames_common <- intersect(colnames, colnames(mat))

  zeromat[rownames_common, colnames_common] <-
    mat[rownames_common, colnames_common] %>%
    Matrix(sparse = TRUE)

  zeromat[rownames, colnames]
}


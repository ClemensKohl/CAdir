
#' Calculate the norm of a row-vector.
#' @param x Row-vector or a matrix of row-vectors.
#' @returns
#' The norm of the vector.
row_norm <- function(x) {
    if (is.matrix(x)) {
        norm <- sqrt(rowSums(x^2))
    } else if (is.null(dim(x))) {
        norm <- sqrt(sum(x^2))
    } else {
        stop("Uknown object.")
    }

    return(norm)
}




#' return k random row indices from a matrix.  
#' @param points Row-wise matrix of points.
#' @param k Number of indices to return.
#' @returns
#' Vector of row-indices.
rand_idx <- function(points, k) {
    idxs <- sample(seq_len(nrow(points)), size = k)
    return(idxs)
}

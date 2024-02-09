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




#' Calculate the cosine of the angle between two vectors.
#' @param a Vector or matrix of row-vectors.
#' @param b Vector or matrix of row-vectors.
#' @returns
#' The cosine of the angle between the vectors.
cosine <- function(a, b) {
    norma <- row_norm(a)
    normb <- row_norm(b)

    if (is.null(dim(a)) && is.null(dim(b))) {
        costheta <- (a %*% b) / (norma * normb)
        costheta <- min(max(costheta, -1), 1)
    } else {
        if (!is.null(dim(b))) b <- t(b)

        costheta <- (a %*% b) / (outer(norma, normb, FUN = "*"))
        costheta <- pmin(pmax(costheta, -1), 1)
    }

    return(costheta)
}

# TODO: Add documentation
get_ang_sim <- function(a, b) {
    sim <- 1 - acos(cosine(a, b)) / pi
    return(sim)
}

# TODO: Add documentation
get_angle <- function(a, b) {
    angle <- acos(cosine(a, b))
    return(angle)
}

# TODO: Add documentation
get_ang_sim <- function(a, b) {
    sim <- 1 - acos(cosine(a, b)) / pi
    return(sim)
}

# TODO: Add documentation
ang_sim_to_deg <- function(ang_sim) {
    deg <- rad2deg((pi - ang_sim * pi))
    return(deg)
}

# TODO: Add documentation
deg2rad <- function(x) {
    r <- (x * pi) / 180
    return(r)
}

# TODO: Add documentation
rad2deg <- function(x) {
    d <- (x * 180) / pi
    return(d)
}

# TODO: Add documentation.
slope <- function(lines, dims = 1:2) {
    if (is.null(nrow(lines)) == 1) {
        d <- lines[dims[2]] / lines[dims[1]]
    } else {
        d <- lines[, dims[2]] / lines[, dims[1]]
    }
    return(d)
}

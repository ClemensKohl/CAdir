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

#' Calculate angular similarity between two vectors.
#' @param a Vector or matrix of row-vectors.
#' @param b Vector or matrix of row-vectors.
#' @returns scalar. The angular similarity.
get_ang_sim <- function(a, b) {
    sim <- 1 - acos(cosine(a, b)) / pi
    return(sim)
}

#' Calculate angle between two vectors.
#' @param a Vector or matrix of row-vectors.
#' @param b Vector or matrix of row-vectors.
#' @returns scalar. The angle between the vectors in radian.
get_angle <- function(a, b) {
    angle <- acos(cosine(a, b))
    return(angle)
}

#' Convert the angular similarity to degrees.
#' Function for manual checking and debugging.
#' @param ang_sim Angular similarity.
#' @returns scalar. The angle in degrees.
ang_sim_to_deg <- function(ang_sim) {
    deg <- rad2deg((pi - ang_sim * pi))
    return(deg)
}

#' Convert the angular similarity to degrees.
#' Function for manual checking and debugging.
#' @param ang_sim Angular similarity.
#' @returns scalar. The angle in degrees.
rad_to_ang_sim <- function(rad) {
    ang_sim <- 1-(rad/pi)
    return(ang_sim)
}


#' Convert degrees to radians.
#' @param x Angle in degrees.
#' @returns Angle in radians.
deg2rad <- function(x) {
    r <- (x * pi) / 180
    return(r)
}

#' Convert radians to degrees.
#' @param x Angle in radians.
#' @returns Angle in degrees.
rad2deg <- function(x) {
    d <- (x * 180) / pi
    return(d)
}

#' Calculate the slope of lines/vectors.
#' Only works in 2D!
#' @param lines Row-wise matrix of lines/vectors.
#' @param dims The 2 dimensions which should be used to calculate the slope.
#' @returns Vector of slopes.
slope <- function(lines, dims = 1:2) {
    if (is.null(nrow(lines)) == 1) {
        d <- lines[dims[2]] / lines[dims[1]]
    } else {
        d <- lines[, dims[2]] / lines[, dims[1]]
    }
    return(d)
}

# ** Plan: **
# There is one main function that performs clustering by direction.
# Depending on the choice of a paremter we branch into Salpha clustering or
# regular splitmerging (e.g. setting cutoff = auto or similar)
#


#' Initialize directions by kmeans++ method
#' @param points Row-wise matrix of points to be clustered.
#' @param k Number of clusters.
#' @returns
#' Vector of indices of the initial directions.
kmeanspp_init <- function(points, k) {
    pnorm <- row_norm(points)
    npoints <- nrow(points)

    center_ids <- integer(k)
    center_ids[1] <- sample.int(npoints, 1)

    for (ii in 2:k) {
        lines <- points[center_ids, ] / pnorm[center_ids]
        ldist <- dist_to_line(points, lines, pnorm)

        probs <- apply(ldist, 1, min)
        probs[center_ids] <- 0

        center_ids[ii] <- sample.int(npoints, 1, prob = probs)
    }

    return(center_ids)
}

#' Calculate the distance of points to the directions
#' @param points Row-wise matrix of points.
#' @param lines Row-wise matrix of lines/directions.
#' @param pnorm Vector of the norm of the points.
#' @returns
#' Matrix of distances of points to the directions.
dist_to_line <- function(points, lines, pnorm) {
    if (is.matrix(lines)) {
        lines <- t(lines)
    }
    proj <- (points %*% lines)
    dist <- pnorm^2 - proj^2
    dist[dist < 0] <- 0
    dist <- sqrt(dist)

    return(dist)
}

# TODO: Add documentation
total_least_squares <- function(points) {

    if (nrow(points) == 1) {
        reg_line <- points / row_norm(points)
    } else {
        suppressWarnings({
            reg_line <- irlba::irlba(points, nv = 1, right_only = TRUE)$v
        })
    }

    return(reg_line)
}

#' Update directions through Total Least Squares Regression
#' @param points Row-wise matrix of points.
#' @param clusters Vector of cluster assignments for the points.
#' @param lines Row-wise matrix of lines.
#' @param k Number of clusters.
#' @returns
#' Matrix of updated directions.
update_line <- function(points, clusters, lines, k) {
    # SVD gives total least squares results
    k <- sort(unique(clusters))
    for (c in seq_len(length(k))) {
        sel <- which(clusters == k[c])

        if (length(sel) == 0){
            next
        } else {
            lines[c, ] <- total_least_squares(points[sel, , drop = FALSE])
        }
    }
    return(lines)
}

#' Basic clustering by directions in CA space
#' @param points Matrix of points to be clustered.
#'  Points are expected to be rows.
#' @param k Number of clusters.
#' @param epochs Number of iterations.
#' @param init Initialization method for the lines. Options are 'rand' and
#' 'kmeanspp'.
#' @param lines Optional. Row-wise matrix of lines.
#'  If not NULL, then the provided lines are used for initialization.
#' @param log Logical. If TRUE, then the function returns a list with
#' the computed distances and directions at each step.
#'
#' @returns
#' An element of class cadir.
#'
#' @export
dirclust <- function(
    points,
    k,
    epochs = 10,
    init = "kmeanspp",
    lines = NULL,
    log = FALSE
) {

    pnorm <- row_norm(points)

    if (is.null(lines)) {
        if (init == "rand") {
            # initialize line
            idxs <- rand_idx(points, k)
            lines <- points[idxs, ]
            lines <- lines / row_norm(lines)
        } else if (init == "kmeanspp") {
            center_ids <- kmeanspp_init(points = points, k = k)
            lines <- points[center_ids, ]
        } else {
            stop("Please use init = 'rand' or 'kmeanspp'.")
        }
    } else {
        # ensure that they are unit vectors.
        lines <- lines / row_norm(lines)
    }

    rownames(lines) <- paste0("line", seq_len(nrow(lines)))

    if (isTRUE(log)) {
        dist_log <- vector(mode = "list", length = epochs)

        dir_log <- vector(mode = "list", length = epochs)
        names(dir_log) <- as.character(seq_len(length(dir_log)))
        dir_log[["0"]] <- lines
    }

    for (i in seq_len(epochs)) {
        # calculate distance to line.
        ldist <- dist_to_line(points, lines, pnorm)

        # find closest line
        clusters <- apply(ldist, 1, which.min)

        # update lines
        lines <- update_line(points, clusters, lines, k)

        if (isTRUE(log)) {
            dir_log[[as.character(i)]] <- lines

            dist_log[[i]] <- vapply(
                seq_len(nrow(ldist)),
                function(x) {
                    ldist[x, clusters[x]]
                },
                numeric(1)
            )
        }
    }

    if (isTRUE(log)) {
        dist_log <- as.data.frame(do.call(rbind, dist_log))
        dimnames(dist_log) <- list(
            paste0("epoch_", seq_len(nrow(dist_log))),
            rownames(points)
        )
        dist_log$epoch <- seq_len(nrow(dist_log))
        log_list <- list("dist_log" = dist_log, "dir_log" = dir_log)
    } else {
        log_list <- list()
    }

    # Remove lines without clusters
    uni_clust <- sort(unique(clusters))
    if (length(uni_clust) != nrow(lines)) {
        lines <- lines[uni_clust, , drop = FALSE]
        ldist <- ldist[, uni_clust, drop = FALSE]
    }

    out <- new("cadir",
        cell_clusters = factor(clusters),
        directions = lines,
        distances = ldist,
        parameters = list(
            "k" = k,
            "epochs" = epochs,
            "init" = init,
            "log" = log
        ),
        log = log_list
    )


    return(out)
}

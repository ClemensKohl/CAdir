#' return k random row indices from a matrix.
#' @param points Row-wise matrix of points.
#' @param k Number of indices to return.
#' @returns
#' Vector of row-indices.
rand_idx <- function(points, k) {
    idxs <- sample(seq_len(nrow(points)), size = k)
    return(idxs)
}

#' Rename clusters and directions so that they match.
#' @param cadir A cadir object.
#' @details
#' rename_clusters takes a cadir object and renames the clusters and directions
#' such that they are again coherent (e.g. from 1:5).
#' @returns
#' A cadir object with renamed clusters and directions.
rename_clusters <- function(cadir) {
    uni_clust <- sort(unique(cadir@clusters))
    nms <- names(cadir@clusters)

    dir_num <- as.numeric(gsub("line", "", rownames(cadir@directions)))

    # Remove directions without any points clustered.
    have_smpls <- which(dir_num %in% uni_clust)
    cadir@directions <- cadir@directions[have_smpls, ]

    # Rename the clusters according to the order of their direction.
    cadir@clusters <- match(cadir@clusters, uni_clust)
    stopifnot(!any(is.na(cadir@clusters)))
    names(cadir@clusters) <- nms

    # Rename the directions according to the cluster numbers.
    rownames(cadir@directions) <- paste0("line", sort(unique(cadir@clusters)))

    return(cadir)
}

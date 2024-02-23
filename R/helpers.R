#' return k random row indices from a matrix.
#' @param points Row-wise matrix of points.
#' @param k Number of indices to return.
#' @returns
#' Vector of row-indices.
rand_idx <- function(points, k) {
    idxs <- sample(seq_len(nrow(points)), size = k)
    return(idxs)
}

#TODO: Fix distances if they are there.

#' Rename clusters and directions so that they match.
#' @param cadir A cadir object.
#' @details
#' rename_clusters takes a cadir object and renames the clusters and directions
#' such that they are again coherent (e.g. from 1:5).
#' @returns
#' A cadir object with renamed clusters and directions.
rename_clusters <- function(cadir) {

    uni_clust <- sort(unique(c(cadir@cell_clusters, cadir@gene_clusters)))
    cell_nms <- names(cadir@cell_clusters)
    gene_nms <- names(cadir@gene_clusters)

    dir_num <- as.numeric(gsub("line", "", rownames(cadir@directions)))

    # Remove directions without any points clustered.
    have_smpls <- which(dir_num %in% uni_clust)
    dir_num <- dir_num[have_smpls]
    cadir@directions <- cadir@directions[have_smpls, , drop = FALSE]

    # Rename the directions according to the cluster numbers.
    new_dir_num <- match(dir_num, uni_clust)
    rownames(cadir@directions) <- paste0("line", new_dir_num)

    if (!is.empty(cadir@distances)) {
        cadir@distances <- cadir@distances[, have_smpls, drop = FALSE]
        colnames(cadir@distances) <- paste0("line", new_dir_num)
    }

    # Rename the clusters according to the order of their direction.
    tmp_cells <- match(cadir@cell_clusters, uni_clust)
    new_lvls <- sort(unique(tmp_cells))


    if (!is.empty(cadir@gene_clusters)) {
        tmp_genes <- match(cadir@gene_clusters, uni_clust)

        new_lvls <- sort(unique(c(tmp_cells, tmp_genes)))

        cadir@gene_clusters <- factor(
            tmp_genes,
            levels = new_lvls
        )
        stopifnot(!any(is.na(cadir@gene_clusters)))
        names(cadir@gene_clusters) <- gene_nms
    }

    cadir@cell_clusters <- factor(
        tmp_cells,
        levels = new_lvls
    )


    stopifnot(!any(is.na(cadir@cell_clusters)))
    names(cadir@cell_clusters) <- cell_nms

    return(cadir)
}

#' Check if a variable is empty (length = 0 but not NULL)
#' @param x variable to check
#' @return TRUE if variable is empty, FALSE otherwise.
is.empty <- function(x) {
    return(isTRUE(length(x) == 0 & !is.null(x)))
}


#' Print cadir object in console.
#' @param object A cadir object
show_cadir <- function(object) {
    stopifnot(is(object, "cadir"))

    ncells <- length(object@cell_clusters)
    ngenes <- length(object@gene_clusters)
    cat("caclust object with", ncells, "cells and", ngenes, "genes.")

    if (!is.empty(CAbiNet::cell_clusters(object)) &&
        !is.empty(CAbiNet::gene_clusters(object))) {

        stopifnot(identical(
            levels(CAbiNet::cell_clusters(object)),
            levels(CAbiNet::gene_clusters(object))
        ))
        cat("\n")
        cat(length(levels(CAbiNet::gene_clusters(object))), "clusters found.")
        cat("\nClustering results:\n\n")
        df <- data.frame(
            "cluster" = levels(CAbiNet::cell_clusters(object)),
            "ncells" = summary(CAbiNet::cell_clusters(object), maxsum = Inf),
            "ngenes" = summary(CAbiNet::gene_clusters(object), maxsum = Inf)
        )

        print(df, row.names = FALSE, right = FALSE)
    } else {
        cat("\nNo biclustering run yet.\n\n")
    }

}

#' @rdname show_cadir
#' @export
setMethod(
    f = "show",
    signature(object = "cadir"),
    function(object) {
        show_cadir(object)
    }
)

#' Convert a factor to a numeric.
#' @details
#' Assumes that the factors are numbers such as "1".
#' @param f A vector of factors.
#' @returns
#' The factors converted to numbers.
f2n <- function(f) {
    f <- as.numeric(as.character(f))
    stopifnot(is.numeric(f))

    return(f)

}

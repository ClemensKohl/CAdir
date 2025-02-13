#' return k random row indices from a matrix.
#' @param points Row-wise matrix of points.
#' @param k Number of indices to return.
#' @returns
#' Vector of row-indices.
rand_idx <- function(points, k) {
    idxs <- sample(seq_len(nrow(points)), size = k)
    return(idxs)
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
    stopifnot(methods::is(object, "cadir"))

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
fc2n <- function(f) {
    n <- as.numeric(as.character(f))
    stopifnot(is.numeric(n))

    return(n)
}

#' Convert factors to numeric.
#' @param f factor vector.
#' @returns
#' A numeric vector.
f2n <- function(f) {
    n <- as.numeric(f)
    stopifnot(is.numeric(n))

    return(n)
}

#' Convert factors to character vector.
#' @param f A vector of factors.
#' @returns
#' A character vector.
f2c <- function(f) {
    c <- as.character(f)
    stopifnot(is.character(c))

    return(c)
}

#' Convert a numeric (or character) to a factor
#' @param n Numeric or character vector.
#' @param lvls The levels for the new factor vector.
x2f <- function(n, lvls) {
    f <- factor(n, levels = lvls)
    return(f)
}



#' Converts a `cadir` object to `Biclust` object
#'
#' @param cadir A cadir object with cell and gene clusters.
#'
#' @return
#' An object of type "Biclust".
#'
#' @export
cadir_to_biclust <- function(cadir) {
    cell_clusters <- cadir@cell_clusters
    gene_clusters <- cadir@gene_clusters
    params <- cadir@parameters

    ctypes <- sort(unique(cell_clusters))
    gtypes <- sort(unique(gene_clusters))
    bitypes <- union(ctypes, gtypes)

    nr <- length(bitypes)

    if (nr == 0) {
        number_x_col <- matrix(0)
        row_x_number <- matrix(0)
    } else {
        number_x_col <- do.call(rbind, lapply(bitypes, function(x) {
            cell_clusters == x
        }))
        row_x_number <- do.call(cbind, lapply(bitypes, function(x) {
            gene_clusters == x
        }))
    }

    rownames(row_x_number) <- names(gene_clusters)
    colnames(row_x_number) <- paste0("BC", bitypes)

    rownames(number_x_col) <- paste0("BC", bitypes)
    colnames(number_x_col) <- names(cell_clusters)

    bic <- methods::new("Biclust",
        "Parameters" = params,
        "RowxNumber" = row_x_number,
        "NumberxCol" = number_x_col,
        "Number" = nr,
        "info" = list("Generated from cell and gene clusters.")
    )

    return(bic)
}

#' Checks if APL S-alpha cutoff is already calculated.
#' @param cadir Cadir object
#' @param fun_args Arguments with which the parent function was called.
#' @returns
#' TRUE if S-alpha store is already calculated. FALSE otherwise.
is_stored <- function(cadir, fun_args) {
    !is.null(cadir@parameters$sa_cutoff) &&
        identical(fun_args$apl_cutoff_reps, cadir@parameters$apl_cutoff_reps) &&
        identical(fun_args$apl_quant, cadir@parameters$call$apl_quant) &&
        identical(fun_args$method, cadir@parameters$call$method)
}

#' Log or append cell clusters to the logs.
#' @param log A (empty) list of previous iterations.
#' @param cadir A CAdir object with cell clusters.
#' @param name Name of the iteration
#' @returns
#' List of iterations with cell clusters.
log_iter <- function(log, cadir, name) {
    log$clusters <- cbind(
        log$clusters,
        stats::setNames(
            data.frame(cadir@cell_clusters),
            name
        )
    )


    log$directions <- rbind(
        log$directions,
        cbind(
            iter = name,
            dirname = rownames(cadir@directions),
            as.data.frame(cadir@directions)
        )
    )

    return(log)
}

#' Convert numeric clusters into a cluster name.
#' @param i Cluster index.
#' @returns
#' Adds prefix "clusters_" to the cluster name.
cl2nm <- function(i) {
    paste0("cluster_", i)
}

#' Searches the dict entry for a given cluster index.
#' @param dict The dictionary (cadir@dict).
#' @param query A string or index corresponding to the name of the directions.
#' @returns
#' The name of the corresponding cluster.
search_dict <- function(dict, query) {
    # names(dict)[dict %in% query]
    names(dict)[base::match(query, dict)]
}

#' Checks if the cluster names confirm to the standard naming.
#' @param nm The cluster name(s) to check.
#' @returns
#' TRUE/FALSE
is_std_name <- function(nm) {
    grepl("cluster_[[:digit:]]+$", nm)
}

#' Extracts the cluster number from standard naming.
#' @param nm A vector of cluster names starting with "cluster_"
#' @returns
#' A vector of cluster numbers.
get_std_num <- function(nm) {
    as.numeric(gsub("^cluster_", "", nm))
}

#' Get indices of cells belonging to a cluster.
#' @param cadir A CAdir object.
#' @param cluster A cluster name.
#' @returns
#' Indices of cells that belong to `cluster`.
get_cluster_idxs <- function(cadir, cluster) {
    stopifnot(is.character(cluster))
    which(cadir@cell_clusters == cluster)
}


#' Checks whether cluster directions have converged.
#' @param now Current CAdir object.
#' @param prev CAdir object from previous iteration.
#' @param cutoff Mean angle between directions from different iterations in
#' degrees below or equal which convergence is reached.
#' @return
#' TRUE or FALSE, depending if convergence criteria is met.
check_convergence <- function(now, prev, cutoff = 0.001) {

    if (nrow(now@directions) != nrow(prev@directions)) {
        return(FALSE)
    }

    ang <- get_angle(now@directions, prev@directions)
    # FIXME: Change from mean to any angle on the diag!
    mean_ang <- mean(diag(ang))
    if (mean_ang <= deg2rad(cutoff)) {
        return(TRUE)
    }

    return(FALSE)
}

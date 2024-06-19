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

# TODO: Add documentation
f2n <- function(f) {
    n <- as.numeric(f)
    stopifnot(is.numeric(n))

    return(n)
}

# TODO: Add documentation.
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

    bic <- new("Biclust",
        "Parameters" = params,
        "RowxNumber" = row_x_number,
        "NumberxCol" = number_x_col,
        "Number" = nr,
        "info" = list("Generated from cell and gene clusters.")
    )

    return(bic)
}

# TODO: Add documentation
#' Checks if APL S-alpha cutoff is already calculated.
#' @param cadir Cadir object
#' @param fun_args Arguments with which the parent function was called.
is_stored <- function(cadir, fun_args) {
    !is.null(cadir@parameters$sa_cutoff) &&
        identical(fun_args$apl_cutoff_reps, cadir@parameters$apl_cutoff_reps) &&
        identical(fun_args$apl_quant, cadir@parameters$call$apl_quant) &&
        identical(fun_args$method, cadir@parameters$call$method)
}

# TODO: Add documentation
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

# TODO: Add Documentation
cl2nm <- function(i) {
    paste0("cluster_", i)
}

# TODO: Add documentation.
search_dict <- function(dict, query) {
    names(dict)[dict %in% query]
}

# TODO: Add documentation.
is_std_name <- function(nm) {
    grepl("cluster_[[:digit:]]+$", nm)
}

# TODO: Add documentation.
get_std_num <- function(nm) {
    as.numeric(gsub("^cluster_", "", nm))
}

#TODO: Add documentaiton
#Possibly delete
get_cluster_idxs <- function(cadir, cluster) {
    stopifnot(is.character(cluster))
    which(cadir@cell_clusters == cluster)
}


# TODO: Find out if you can have the class without cabinet package.


#' Checks if cadir-class was constructed correctly
#'
#' @param object cadir object
#'
#' @return
#' If object is a valid cadir object returns TRUE, otherwise the errors.
#'
check_cadir <- function(object) {
    stopifnot(is(object, "cadir"))

    errors <- character()

    ndir <- nrow(object@directions)
    n_cell_cl <- unique(object@cell_clusters)
    n_gene_cl <- unique(object@gene_clusters)

    n_cl <- length(unique(n_cell_cl, n_gene_cl))

    if (ndir != n_cl) {
        msg <- c(
            "Number of cell & gene clusters not equal to number of directions."
        )
        errors <- c(errors, msg)
    }

    if (ndir != ncol(object@distances) || ndir != ncol(object@distances)) {
        msg <- c("Distances dont match directions.")
        errors <- c(errors, msg)
    }

    if (is.null(names(object@cell_clusters))) {
        msg <- c("Cell clusters have no names.")
        errors <- c(errors, msg)
    }

    if (!is.empty(object@gene_clusters) && is.null(names(object@gene_clusters))) {
        msg <- c("Gene clusters have no names.")
        errors <- c(errors, msg)
    }

    if (length(errors) == 0) TRUE else errors
}


#' An S4 class for the CA directional clustering.
#' @name cadir-class
#' @rdname cadir-class
#' @description
#' Class to store biclustering by directions results.
#' @slot cell_clusters factors. The assigned cell clusters with cell names in
#' the names attribute.
#' @slot gene_clusters factors. The assigned gene clusters with gene names in
#' the names attribute.
#' @slot directions Matrix of directions by which the data was clustered.
#' @slot distances Matrix of distances of points to the respective directions.
#' @slot parameters List of used parameters and function name with which results
#' were generated.
#' @export
setClass("cadir",
    contains = "caclust",
    representation(
        cell_clusters = "factor",
        gene_clusters = "factor",
        directions = "matrix",
        distances = "matrix",
        log = "list",
        parameters = "list",
        plots = "list"
    ),
    prototype(
        cell_clusters = factor(),
        gene_clusters = factor(),
        directions = matrix(0, 0, 0),
        distances = matrix(0, 0, 0),
        parameters = list(),
        log = list(),
        plots = list("splits" = list(),
                     "merges" = list(),
                     "clusters" = list())
    ),
    validity = check_cadir
)



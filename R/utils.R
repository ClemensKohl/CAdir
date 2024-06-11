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
    #FIXME: Likely fails with annotated dataset.
    #FIXME: Change line to cl2dir mechanism or cluster_

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
f2n <- function(f) {
    n <- as.numeric(as.character(f))
    stopifnot(is.numeric(n))

    return(n)
}

# TODO: Add documentation.
f2c <- function(f) {
    c <- as.character(f)
    stopifnot(is.character(c))

    return(c)
}

# TODO: Add documentation
#' Convert a numeric (or character) to a factor
n2f <- function(n, lvls) {
    f <- factor(n, levels = lvls)
    return(f)
}

#' Build a sub-graph that can be joined to a larger graph.
#' By itself it is a valid graph too.
#' @param before A vector (optimally of factors) with all cluster/node assignments.
#' @param after A vector of same length and order as `before`,
#' with cluster/node assignments after some change.
#' @param before_nm The name of the node before the change.
#' @param after_nm The name of the node after the change.
#' @returns
#' A data.frame with two columns ("from", "to") representing the edges in the graph.
#' The data frame can be used to generate a graph.
build_sub_graph <- function(before,
                            after,
                            before_nm = "start",
                            after_nm = "end") {
    graph <- data.frame()
    node <- unique(before)

    for (n in node) {
        n_nm <- paste0(before_nm, "-", n)
        cluster <- which(before == n)
        new_clusters <- unique(after[cluster])

        for (nn in new_clusters) {
            nn_nm <- paste0(after_nm, "-", nn)
            graph <- rbind(graph, data.frame(from = n_nm, to = nn_nm))
        }
    }
    return(graph)
}

#' Build igraph from a cadir results.
#' @param cadir A cadir object with valid clustering results.
#' @param rm_redund If TRUE, removes all clustering iterations where nothing changed (no splits/merges).
#' @returns
#' A directed igraph object with all splits and merges.
build_graph <- function(cadir, rm_redund = FALSE, keep_end = TRUE) {
    graph_list <- list()
    cls <- cadir@log$clusters

    cls_nodes <- colnames(cls)
    node_pattern <- c("iter_0|split|merge|end")
    sel <- which(grepl(node_pattern, cls_nodes))

    # We loop through all the nodes that contain the pattern in their name and
    # compare the new clusters to the ones from the intermediate dirclust calls.
    # This way we only record in the graph what actually changes during the
    # split/merge,
    # not the arbitrary shifting during the cluster refinement. However,
    # we need to change the name of the previous graph node to that of the
    # last split/merge in order to make
    # a coherent graph.
    for (i in seq_len(length(sel))) {
        bef_cls <- cls[, sel[i] - 1]
        aft_cls <- cls[, sel[i]]

        # FIXME: There is a bug with this.
        is_end <- (i == length(sel) && keep_end)
        if (!is_end &&
            all(bef_cls == aft_cls) &&
            isTRUE(rm_redund)) next

        graph_list[[i]] <- build_sub_graph(
            before = bef_cls,
            after = aft_cls,
            before_nm = ifelse(i == 1, "root", cls_nodes[sel[last_i]]),
            after_nm = cls_nodes[sel[i]]
        )
        last_i <- i
    }

    graph <- do.call(rbind, graph_list)
    rownames(graph) <- NULL
    graph <- igraph::graph_from_data_frame(graph, directed = TRUE)

    return(graph)
}

# FIXME: Prevent opening plot when calling it! ggplotGrob

#' Get relative x and y values in relation to the plotting panel from a ggplot.
#' Solution adapted from Anwer by Allan Cameron at:
#' https://stackoverflow.com/a/60857307/1376616
#' @param gg_plot A ggplot object.
#' @returns
#' dataframe with x and y values relative to the plotting panel.
get_x_y_values <- function(gg_plot) {
    img_dim <- grDevices::dev.size("cm") * 10
    gt <- ggplot2::ggplotGrob(gg_plot)
    to_mm <- function(x) grid::convertUnit(x, "mm", valueOnly = TRUE)
    n_panel <- which(gt$layout$name == "panel")
    panel_pos <- gt$layout[n_panel, ]
    panel_kids <- gtable::gtable_filter(gt, "panel")$grobs[[1]]$children
    point_grobs <- panel_kids[[grep("point", names(panel_kids))]]
    from_top <- sum(to_mm(gt$heights[seq(panel_pos$t - 1)]))
    from_left <- sum(to_mm(gt$widths[seq(panel_pos$l - 1)]))
    from_right <- sum(to_mm(gt$widths[-seq(panel_pos$l)]))
    from_bottom <- sum(to_mm(gt$heights[-seq(panel_pos$t)]))
    panel_height <- img_dim[2] - from_top - from_bottom
    panel_width <- img_dim[1] - from_left - from_right
    xvals <- as.numeric(point_grobs$x)
    yvals <- as.numeric(point_grobs$y)
    yvals <- yvals * panel_height + from_bottom
    xvals <- xvals * panel_width + from_left
    data.frame(x = xvals / img_dim[1], y = yvals / img_dim[2])
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

log_iter <- function(log, cadir, name) {

    log$clusters <- cbind(
        log$clusters,
        stats::setNames(
            data.frame(f2n(cadir@cell_clusters)),
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

cl2nm <- function(i) {
    paste0("cluster_", i)
}

set_clusters <- function(clusters, directions, dict) {

}

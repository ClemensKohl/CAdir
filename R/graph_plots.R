# FIXME: Prevent opening plot when calling it! ggplotGrob
#
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


#' Build a sub-graph that can be joined to a larger graph.
#' By itself it is a valid graph too.
#' @param before A vector (optimally of factors) with
#' all cluster/node assignments.
#' @param after A vector of same length and order as `before`,
#' with cluster/node assignments after some change.
#' @param before_nm The name of the node before the change.
#' @param after_nm The name of the node after the change.
#' @returns
#' A data.frame with two columns ("from", "to") representing
#' the edges in the graph.
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
#' @param rm_redund If TRUE, removes all clustering iterations
#' where nothing changed (no splits/merges).
#' @param keep_end Adds the last iteration, even if no change to previous.
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

        is_end <- (i == length(sel) && keep_end)
        if (!is_end &&
            all(f2c(bef_cls) == f2c(aft_cls)) &&
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


#' Makes a simple plot of the splits and merges.
#' @param cadir A cadir object with valid clustering results.
#' @param rm_redund If TRUE, only shows an iteration if something changes.
#' @param alpha Between 0 and 1. Sets opacity of the nodes.
#' @param size Size of the nodes.
#' @returns A plot that that shows the
#' splits and merges as decisions in a graph.
plot_sm_graph <- function(cadir,
                          rm_redund = TRUE,
                          size = 3,
                          alpha = 1) {
    graph <- build_graph(cadir = cadir, rm_redund = rm_redund)

    lgraph <- ggraph::create_layout(graph, layout = "tree")

    ggraph::set_graph_style(plot_margin = ggplot2::margin(0, 0, 0, 0))
    p <- ggraph::ggraph(lgraph) +
        ggraph::geom_edge_link(color = mpi_extend_pal()(2)[2]) +
        ggraph::geom_node_point(size = size, alpha = alpha, color = mpi_extend_pal()(1)) +
        ggraph::theme_graph()

    return(p)
}

#' Plots the graph of the clustering splits and merges
#' and overlays APL plots over the nodes.
#' @inheritParams cluster_apl
#' @inheritParams plot_sm_graph
#' @inheritParams build_graph
#' @param caobj A cacomp object.
#' @param org Organisme to use for the automatic annotation.
#' Either "mm" for mouse or "hs" for human.
#' @param annotate_clusters If TRUE uses automatic cell type annotation via
#' `annotate_biclustering` to annotate the clusters in each level.
#' @returns
#' A ggplot object showing the split-merge graph and APL plots for each cluster.
#' @export
sm_plot <- function(cadir,
                    caobj,
                    rm_redund = TRUE,
                    show_cells = TRUE,
                    show_genes = FALSE,
                    highlight_cluster = FALSE,
                    annotate_clusters = FALSE,
                    org = "mm",
                    keep_end = TRUE,
                    inlet_side = 0.08) {
    # TODO: Simplify function.
    base::stopifnot(
        "Set either `show_cells` or `show_genes` to TRUE." =
            isTRUE(show_cells) || isTRUE(show_genes)
    )

    graph <- build_graph(
        cadir = cadir,
        rm_redund = rm_redund,
        keep_end = keep_end
    )

    lgraph <- ggraph::create_layout(graph, layout = "tree")

    ggraph::set_graph_style(plot_margin = ggplot2::margin(0, 0, 0, 0))
    bg <- ggraph::ggraph(lgraph) +
        ggraph::geom_edge_link() +
        ggraph::geom_node_point(alpha = 1)

    bg_coords <- get_x_y_values(bg)

    cls <- cadir@log$clusters
    dirs <- cadir@log$directions

    nodes <- names(igraph::V(graph))

    old_iter_nm <- ""
    for (i in seq_len(nrow(lgraph))) {
        node_nm <- nodes[i]

        name_elems <- base::strsplit(node_nm, "-", fixed = TRUE)[[1]]
        # name_elems <- stringr::str_split_1(node_nm, "-")

        if (name_elems[1] == "root") next

        iter_nm <- name_elems[1]
        cluster <- name_elems[2]

        grp_idx <- base::which(cls[, iter_nm] == cluster)

        is_iter_dirs <- dirs$iter == iter_nm
        coord_column <- !colnames(dirs) %in% c("iter", "dirname")

        tmp_dirs <- dirs[is_iter_dirs, coord_column]
        rownames(tmp_dirs) <- dirs[is_iter_dirs, "dirname"]

        cluster_idx <- base::which(rownames(tmp_dirs) == cluster)
        dir <- tmp_dirs[cluster_idx, ]

        if (iter_nm != old_iter_nm) {
            tmp_ccs <- x2f(cls[, iter_nm])
            names(tmp_ccs) <- rownames(caobj@prin_coords_cols)

            tmp_cadir <- methods::new(
                "cadir",
                cell_clusters = tmp_ccs,
                directions = as.matrix(tmp_dirs)
            )

            if (is.null(cadir@parameters$qcutoff)) {
                cadir@parameters$qcutoff <- 0.8
            }

            tmp_cadir@gene_clusters <- assign_genes(
                caobj = caobj,
                cadir = tmp_cadir,
                qcutoff = cadir@parameters$qcutoff
            )

            if (isTRUE(annotate_clusters)) {
                suppressWarnings({
                    tmp_cadir <- annotate_biclustering(
                        obj = tmp_cadir,
                        universe = rownames(caobj@std_coords_rows),
                        org = org,
                        alpha = 0.05,
                        min_size = 10,
                        max_size = 500
                    )
                })
            }
            old_iter_nm <- iter_nm
        }

        cluster <- rownames(tmp_cadir@directions)[cluster_idx]
        rownames(dir) <- cluster

        # colour_by_group <- !highlight_cluster

        p <- cluster_apl(
            caobj = caobj,
            cadir = tmp_cadir,
            direction = as.numeric(dir),
            group = grp_idx,
            cluster = cluster,
            show_cells = show_cells,
            show_genes = show_genes,
            highlight_cluster = highlight_cluster,
            show_lines = FALSE,
            point_size = 0.3
        )
        if (isTRUE(annotate_clusters)) {
            p <- p +
                ggplot2::ggtitle(cluster) +
                theme_blank(
                    title = ggplot2::element_text(color = "black",
                                                  size = 10, face = "bold"),
                    text = ggplot2::element_text()
                )
        } else {
            #TODO: We need to pick a color palette for a large number of clusters
            # scale_color_mpimg(name = "mpimg") +
            p <- p + theme_blank()
        }

        # Ensure that we dont plot outside of the window:
        lradj <- tbadj <- 0
        mid_length <- inlet_side/2
        if (bg_coords[i, 1] < mid_length) {
            lradj <- mid_length - bg_coords[i, 1]
        }
        if ((bg_coords[i, 1] + mid_length) > 1) {
            lradj <- 1 - (mid_length + bg_coords[i, 1])
        }
        if (bg_coords[i, 2] < mid_length) {
            tbadj <- mid_length - bg_coords[i, 2]
        }
        if ((bg_coords[i, 2] + mid_length) > 1) {
            tbadj <- 1 - (mid_length + bg_coords[i, 2])
        }


        bg <- bg +
            patchwork::inset_element(p,
                left = bg_coords[i, 1] - mid_length + lradj,
                right = bg_coords[i, 1] + mid_length + lradj,
                top = bg_coords[i, 2] + mid_length + tbadj,
                bottom = bg_coords[i, 2] - mid_length + tbadj,
                align_to = "panel"
            )
    }

    return(bg)
}

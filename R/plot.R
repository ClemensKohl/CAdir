#' Create a model to project points into an Association Plot
#' @param caobj A cacomp object.
#' @param direction Normed direction vector of the APL plot.
#' @param group A vector of indices which indicate the points
#' that belong to the cluster. Only needed here to orient the plot.
#' @returns
#' A model that can be used to project new points onto the APL plot.
apl_model <- function(
    caobj,
    direction,
    group = NULL) {

    stopifnot(methods::is(caobj, "cacomp"))

    cent <- caobj@prin_coords_cols

    avg_group_coords <- direction
    length_vector_group <- sqrt(drop(avg_group_coords %*% avg_group_coords))

    # The line sometimes point into the "wrong" direction.
    # We can determine the cosine, if its negative we flip the x coords.

    cosangle <- 1

    if (!is.null(group)) {
        subgroup <- cent[group, ]

        if (length(group) == 1) {
            group_mean <- subgroup # single sample
        } else {
            group_mean <- colMeans(subgroup) # centroid vector.
        }

        cosangle <- cosine(a = group_mean, b = direction)

        # group_norm <- row_norm(group_mean)
        # gx <- drop(group_mean %*% avg_group_coords) / group_norm
        #
        # if (sign(gx) == -1) {
        #     avg_group_coords <- -avg_group_coords
        # }

        if (cosangle < 0) {
            avg_group_coords <- -avg_group_coords
        }
    }


    model <- function(vec) {
        length_vector <- row_norm(vec)
        cordx <- drop(vec %*% avg_group_coords) / length_vector_group
        cordy <- suppressWarnings(sqrt(length_vector^2 - cordx^2))

        cordx[is.na(cordx)] <- 0
        cordy[is.na(cordy)] <- 0

        return(cbind("x" = cordx, "y" = cordy))
    }

    return(model)
}

#' Plot a cluster with the respective direction/line in an APL.
#' @inheritParams apl_model
#' @param cadir A cadir object for which to compute the APL
#' @param direction Direction of the APL plot.
#' @param cluster The cluster (if any) to highlight cells and genes by.
#' @param group Group to determine the correct direction of the APL plot.
#' Can be the same as the cluster.
#' @param show_cells If TRUE, points (cells) are plotted.
#' @param show_genes If TRUE, genes are plotted.
#' @param show_lines If TRUE, the directions in cadir are plotted.
#' @param highlight_cluster If TRUE, highlights the points in belonging to
#'  `cluster`, unless colour_by_group is TRUE.
#' @param colour_by_group If TRUE highlights cells that belong to `group`.
#' @param point_size Size of the points (cells).
#' @returns
#' An APL plot (ggplot2 object).
#' @export
cluster_apl <- function(caobj,
                        cadir,
                        cluster = NULL,
                        direction,
                        group,
                        show_cells = TRUE,
                        show_genes = FALSE,
                        show_lines = TRUE,
                        highlight_cluster = TRUE,
                        colour_by_group = FALSE,
                        point_size = 1.5) {

    # TODO: Check if you can simplify the cluster/group assignments.
    stopifnot(methods::is(caobj, "cacomp"))
    stopifnot(methods::is(cadir, "cadir"))

    if (!is.null(cluster)) cluster <- as.character(cluster)

    all_cls <- as.character(
        sort(unique(c(cadir@cell_clusters, cadir@gene_clusters)))
    )

    if (!is.null(cluster) && (!cluster %in% all_cls)) cluster <- NULL
    if (is.null(cluster)) {
        # Kinda redundant. placeholder if I want to do deal with special case.
        ccluster <- NULL
        gcluster <- NULL

        cell_grp <- seq_len(length(cadir@cell_clusters))
        gene_grp <- seq_len(length(cadir@gene_clusters))
    } else {
        # ccluster <- n2f(cluster, lvls = levels(cadir@cell_clusters))
        # gcluster <- n2f(cluster, lvls = levels(cadir@gene_clusters))
        cell_grp <- which(f2c(cadir@cell_clusters) == cluster)
        gene_grp <- which(f2c(cadir@gene_clusters) == cluster)
    }

    # cat("\nCluster:", cluster)
    # cat("\nType of cluster", class(cluster))
    # cat("\n cell clusters:", unique(cadir@cell_clusters))
    # cat("\n cell clusters class:", class(cadir@cell_clusters))
    # cat("\n directions:", nrow(cadir@directions))

    # sel <- match(cluster, sort(unique(cadir@cell_clusters)))
    # direction <- cadir@directions[sel, ]

    # ensure that clusters and directions are coherent
    # cadir <- rename_clusters(cadir)
    if (nrow(cadir@directions) == 2) {
        ang <- min(
            rad2deg(get_angle(cadir@directions[1, ], cadir@directions[2, ])),
            rad2deg(get_angle(-cadir@directions[1, ], cadir@directions[2, ]))
        )
    } else {
        ang <- 0
    }

    model <- apl_model(
        caobj = caobj,
        direction = direction,
        group = group
    )

    dapl <- model(cadir@directions)

    if (show_cells && !show_genes) {
        capl <- model(caobj@prin_coords_cols)
        df <- as.data.frame(capl)
        df$sample <- rownames(df)
        df$type <- "cell"
    } else if (show_cells && show_genes) {
        capl <- model(caobj@std_coords_cols)
        gapl <- model(caobj@prin_coords_rows)
        df <- as.data.frame(capl)
        df$sample <- rownames(df)
        df$type <- "cell"

        dfg <- as.data.frame(gapl)
        dfg$sample <- rownames(dfg)
        dfg$type <- "gene"

        df <- rbind(df, dfg)
    } else if (!show_cells && show_genes) {
        gapl <- model(caobj@prin_coords_rows)
        df <- as.data.frame(gapl)
        df$sample <- rownames(df)
        df$type <- "gene"
    } else {
        df <- data.frame(
            "sample" = c(),
            "type" = c()
        )
    }

    # If the line points into the opposite direction of points
    # we flip the line.
    for (d in seq_len(nrow(dapl))) {
        sel <- match(
            names(cadir@cell_clusters)[as.numeric(cadir@cell_clusters) == d],
            rownames(caobj@prin_coords_cols)
        )

        if (length(sel) == 0) next

        if (length(sel) > 1) {
            grp_mean <- colMeans(model(caobj@prin_coords_cols)[sel, ])
        } else {
            grp_mean <- model(caobj@prin_coords_cols)[sel, ]
        }

        if (sign(grp_mean[1]) != sign(dapl[d, 1])) {
            dapl[d, ] <- c(-1, 1) * dapl[d, ]
        }
    }

    if (isTRUE(highlight_cluster) && (show_cells || show_genes)) {
        if (isTRUE(colour_by_group)) {
            sel <- group
        } else {
            sel <- c()
            if (show_cells) {
                sel_cells <- match(
                    rownames(caobj@prin_coords_cols)[cell_grp],
                    df$sample
                )
                sel <- c(sel, sel_cells)
            }
            if (show_genes) {
                sel_genes <- match(
                    rownames(caobj@prin_coords_rows)[gene_grp],
                    df$sample
                )
                sel <- c(sel, sel_genes)
            }
        }
        sel <- na.omit(sel)

        df$cluster <- "other"
        df$cluster[sel] <- "cluster"
        df$cluster <- factor(df$cluster, levels = c("other", "cluster"))
        ord <- order(df$cluster)
        df <- df[ord, ]
    } else if (show_cells || show_genes){
        df$cluster <- 0

        sel_cells <- match(names(cadir@cell_clusters), df$sample)
        sel_genes <- match(names(cadir@gene_clusters), df$sample)
        sel_cells <- na.omit(sel_cells)
        sel_genes <- na.omit(sel_genes)

        df$cluster[sel_cells] <- cadir@cell_clusters
        df$cluster[sel_genes] <- cadir@gene_clusters

        df$cluster <- factor(df$cluster,
            levels = sort(unique(c(
                0,
                cadir@cell_clusters,
                cadir@gene_clusters
            )))
        )
    }

    # if conditions, conc. type and cluster.
    if (highlight_cluster && show_cells && show_genes) {
        df$cluster <- paste0(df$type, "_", df$cluster)
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = cluster))

    if (isTRUE(show_cells) || isTRUE(show_genes)) {
        if (show_cells && show_genes) {
            size_factor <- 1.5
        } else {
            size_factor <- 1
        }
        p <- p +
            ggplot2::geom_point(ggplot2::aes(shape = type, size = type)) +
            ggplot2::scale_size_manual(values = c(
                "cell" = point_size,
                "gene" = point_size / size_factor
            )) +
            ggplot2::scale_shape_manual(values = c(
                "cell" = 19,
                "gene" = 1
            ))

        if (isTRUE(highlight_cluster)) {
            if (show_cells && show_genes) {
                p <- p + ggplot2::scale_color_manual(values = c(
                    "cell_cluster" = "#c6d325",
                    "gene_cluster" = "#ef7c00",
                    "cell_other" = "#006c66",
                    "gene_other" = "#777777"
                ))
            } else {
                p <- p + ggplot2::scale_color_manual(values = c(
                    "cluster" = "#c6d325",
                    "other" = "#006c66"
                ))
            }
        }
    }


    if (isTRUE(show_lines)) {
        for (d in seq_len(nrow(dapl))) {
            if (isTRUE(all.equal(dapl[d, ],
                c(1, 0),
                tolerance = 1e-4,
                check.attributes = FALSE
            )) ||
                isTRUE(all.equal(dapl[d, ],
                    c(-1, 0),
                    tolerance = 1e-4,
                    check.attributes = FALSE
                ))) {
                lcolor <- "black"
                ltype <- "solid"
            } else {
                lcolor <- "red"
                ltype <- "dashed"

                if (isTRUE(highlight_cluster)) {
                    lcolor <- "#006c66"
                }
            }

            p <- p + ggplot2::geom_abline(
                intercept = 0,
                slope = slope(lines = dapl[d, ], dims = 1:2),
                color = lcolor,
                linetype = ltype,
                size = 1
            ) +
                ggplot2::geom_point(
                    data = data.frame(x = 0, y = 0),
                    ggplot2::aes(x, y), color = lcolor
                )
        }
    }

    p <- p + ggplot2::ggtitle(paste0(
        "Cluster: ",
        as.character(cluster),
        ", CA-angle: ",
        round(ang, 2)
    )) +
        ggplot2::theme_bw()

    return(p)
}

#' Plot for the clustering results wich shows the
#' relationship between the clusters. In the diagonal
#' an APL plot for the respective cluster is shown.
#' @param cadir A cadir object with valid cell clustering results.
#' @param caobj A cacomp object.
#' @inheritParams cluster_apl
#' @returns A plot that summarizes the cell clustering results and
#' how the clusters relate to each other.
#' @export
plot_results <- function(cadir,
                         caobj,
                         highlight_cluster = TRUE,
                         show_cells = TRUE,
                         show_genes = FALSE) {
    size <- 1
    pls <- list()
    cls <- sort(unique(cadir@cell_clusters))

    for (i in seq_along(cls)) {
        for (j in seq_along(cls)) {
            sel <- which(cadir@cell_clusters == cls[i] |
                         cadir@cell_clusters == cls[j])
            sel_dir <- unique(c(f2n(cls[i]), f2n(cls[j])))

            sub_cak <- methods::new("cadir",
                cell_clusters = cadir@cell_clusters[sel],
                directions = cadir@directions[sel_dir, , drop = FALSE]
            )
            # cat("\n\ni,j", i, ",", j)
            if (show_cells) sc <- i == j
            if (show_genes) sg <- i == j
            p <- cluster_apl(
                caobj = caobj,
                cadir = sub_cak,
                direction = cadir@directions[f2n(cls[i]), ],
                group = which(cadir@cell_clusters == cls[i]),
                cluster = cls[i],
                show_cells = sc,
                show_genes = sg,
                show_lines = i != j,
                highlight_cluster = highlight_cluster,
                colour_by_group = TRUE
            ) +
                ggplot2::ggtitle("") +
                ggplot2::theme(
                    legend.position = "none",
                    axis.title.x = ggplot2::element_blank(),
                    axis.text.x = ggplot2::element_blank(),
                    axis.ticks.x = ggplot2::element_blank(),
                    axis.title.y = ggplot2::element_blank(),
                    axis.text.y = ggplot2::element_blank(),
                    axis.ticks.y = ggplot2::element_blank()
                )


            if (i == j) {
                p$layers[[1]]$aes_params$size <- size
            } else {
                p$layers[[2]]$aes_params$size <- size
            }
            pls[[paste0(i, "_", j)]] <- p
        }
    }

    fig <- ggpubr::ggarrange(plotlist = pls,
                             nrow = length(cls),
                             ncol = length(cls))
    return(fig)
}

#' Summarizes the cell clustering results in a single plot.
#' @param cadir A cadir object with valid cell clustering results.
#' @param caobj A cacomp object.
#' @param point_size Size of the points (cells).
#' @inheritParams cluster_apl
#' @returns A plot that summarizes the cell clustering results.
#' @export
plot_clusters <- function(cadir, caobj, point_size = 1, show_genes = FALSE) {
    pls <- list()
    cls <- sort(unique(cadir@cell_clusters))

    for (i in seq_along(cls)) {
        p <- cluster_apl(
            caobj = caobj,
            cadir = cadir,
            direction = cadir@directions[f2n(cls[i]), ],
            cluster = cls[i],
            group = which(cadir@cell_clusters == cls[i]),
            show_cells = TRUE,
            show_genes = show_genes,
            show_lines = FALSE,
            highlight_cluster = TRUE,
            colour_by_group = FALSE
        ) +
            ggplot2::ggtitle(paste0("cluster_", i)) +
            ggplot2::theme(
                legend.position = "none",
                axis.title.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank()
            )

        p$layers[[1]]$aes_params$size <- point_size

        pls[[i]] <- p
    }

    fig <- ggpubr::ggarrange(plotlist = pls, nrow = ceiling(sqrt(length(cls))), ncol = ceiling(sqrt(length(cls))))
    return(fig)
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
        ggraph::geom_edge_link(color = mpi_pal()(2)[2]) +
        ggraph::geom_node_point(size = size, alpha = alpha, color = mpi_pal()(1)) +
        ggraph::theme_graph()

    return(p)
}



#' Creates a sankey plot of the clustering results with all splits and merges.
#' @param cadir A cadir object with valid clustering results.
#' @param rm_redund If TRUE, only shows an iteration if something changes.
#' @returns A sankey plot of the clustering results.
plot_ggsankey <- function(cadir, rm_redund = TRUE) {
    node_pattern <- c("root|iter_0|split|merge|end")
    sel <- which(grepl(node_pattern, colnames(cadir@log$clusters)))
    sub_cls <- cadir@log$clusters[, sel]

    if (isTRUE(rm_redund)) {
        del_cl <- c()
        for (c in seq_len(ncol(sub_cls))) {
            if (c == 1) next
            if (all(sub_cls[, c] == sub_cls[, c - 1])) {
                del_cl <- c(del_cl, c)
            }
        }
        if (length(del_cl) > 0) {
            sub_cls <- sub_cls[, -del_cl]
        }
    }

    sank <- ggsankey::make_long(
        sub_cls,
        colnames(sub_cls)
    )
    p <- ggplot2::ggplot(
        sank,
        ggplot2::aes(
            x = x,
            next_x = next_x,
            node = node,
            next_node = next_node,
            fill = factor(node),
            label = node
        )
    ) +
        ggsankey::geom_sankey(node_color = 1, flow_alpha = 0.7) +
        ggsankey::geom_sankey_label(size = 3.5, color = 1, fill = "white") +
        ggplot2::scale_fill_viridis_d(option = "A", alpha = 0.95) +
        ggsankey::theme_sankey(base_size = 12) +
        ggplot2::theme(legend.position = "none")
    return(p)
}

#' Creates a sankey plot of the clustering results with all splits and merges.
#' @param cadir A cadir object with valid clustering results.
#' @param rm_redund If TRUE, only shows an iteration if something changes.
#' @returns A sankey plot of the clustering results.
#' @export
plot_sankey <- function(cadir, rm_redund = rm_redund) {
    # graph <- build_graph(cadir = cadir, rm_redund = rm_redund)
    #
    # mbms <- igraph::membership(graph)
    # ntwd3 <- networkD3::igraph_to_networkD3(graph)

    node_pattern <- c("root|iter_0|split|merge|end")
    sel <- which(grepl(node_pattern, colnames(cadir@log$clusters)))
    sub_cls <- cadir@log$clusters[, sel]

    if (isTRUE(rm_redund)) {
        del_cl <- c()
        for (c in seq_len(ncol(sub_cls))) {
            if (c == 1) next
            if (all(sub_cls[, c] == sub_cls[, c - 1])) {
                del_cl <- c(del_cl, c)
            }
        }
        if (!is.null(del_cl)) {
            sub_cls <- sub_cls[, -del_cl]
        }
    }

    nodes_nms <- c()
    links <- data.frame()

    stepnm <- colnames(sub_cls)

    for (n in seq_len(ncol(sub_cls))) {
        nodes_nms <- c(nodes_nms, paste0(stepnm[n], "_", unique(sub_cls[, n])))

        if (n == ncol(sub_cls)) next

        from <- unique(sub_cls[, n])

        for (l in seq_len(length(from))) {
            idx <- which(sub_cls[, n] == from[l])
            to <- unique(sub_cls[idx, n + 1])

            for (t in seq_len(length(to))) {
                val <- sum(sub_cls[idx, n + 1] == to[t])

                tmpdf <- data.frame(
                    "source" = paste0(stepnm[n], "_", from[l]),
                    "target" = paste0(stepnm[n + 1], "_", to[t]),
                    "value" = val
                )

                links <- rbind(links, tmpdf)
            }
        }
    }

    nodes_nms <- factor(nodes_nms, levels = nodes_nms)
    node_ids <- as.numeric(nodes_nms) - 1 # 0-idxed
    nodes <- data.frame("id" = node_ids, "name" = nodes_nms)

    links$source <- factor(links$source, levels = nodes_nms)
    links$target <- factor(links$target, levels = nodes_nms)
    links$source_id <- as.numeric(links$source) - 1
    links$target_id <- as.numeric(links$target) - 1


    networkD3::sankeyNetwork(
        Links = links,
        Nodes = nodes,
        Source = "source_id",
        Target = "target_id",
        Value = "value",
        NodeID = "name",
        fontSize = 12,
        nodeWidth = 30,
        NodeGroup = "name",
        # LinkGroup = "source"
    )
}

#' Plots the graph of the clustering splits and merges
#' and overlays APL plots over the nodes.
#' @inheritParams plot_sm_graph
#' @inheritParams cluster_apl
#' @param caobj A cacomp object.
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
                    org = "mm") {

    #FIXME: Genes are basically impossible to tell from cells
    graph <- build_graph(cadir = cadir, rm_redund = rm_redund)

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
        cluster <- as.numeric(name_elems[2])

        grp_idx <- which(cls[, iter_nm] == cluster)

        # TODO: We calculate the directions new. Shouldn't we save them somehow?
        # dir <- total_least_squares(caobj@prin_coords_cols[grp_idx, ])

        is_iter_dirs <- dirs$iter == iter_nm
        dir <- dirs[is_iter_dirs, colnames(dirs) != "iter"]
        # FIXME: this is not reliable. Maybe add the cluster to the logging info.
        dir <- dir[cluster, ]

        if (isTRUE(annotate_clusters)) {
            if (iter_nm != old_iter_nm) {
                tmp_ccs <- n2f(cls[, iter_nm])
                names(tmp_ccs) <- rownames(caobj@prin_coords_cols)

                tmp_cadir <- methods::new(
                    "cadir",
                    cell_clusters = tmp_ccs,
                    directions = as.matrix(dirs[is_iter_dirs,
                                           colnames(dirs) != "iter"])
                )

                tmp_cadir@gene_clusters <- CAdir:::assign_genes(
                    caobj = caobj,
                    cadir = tmp_cadir,
                    qcutoff = 0.8
                )

                tmp_cadir <- CAdir::annotate_biclustering(
                    obj = tmp_cadir,
                    universe = rownames(caobj@std_coords_rows),
                    org = org,
                    alpha = 0.05,
                    min_size = 10,
                    max_size = 500
                )
                old_iter_nm <- iter_nm
            }

            cell_type <- rownames(tmp_cadir@directions)[cluster]
        }

        p <- cluster_apl(
            caobj = caobj,
            cadir = cadir,
            direction = as.numeric(dir),
            group = grp_idx,
            cluster = cluster,
            show_cells = show_cells,
            show_genes = show_genes,
            highlight_cluster = highlight_cluster,
            colour_by_group = TRUE,
            show_lines = FALSE,
            point_size = 0.3
        )
        if (isTRUE(annotate_clusters)) {
            p <- p + ggplot2::ggtitle(cell_type) +
                theme_blank(
                    title = ggplot2::element_text(color = "black", size = 10, face = "bold"),
                    text = ggplot2::element_text()
                )
        } else {
            # scale_color_mpimg(name = "mpimg") + #FIXME: We need to pick a color palette for a large number of clusters
            p <- p + theme_blank()
        }

        bg <- bg +
            patchwork::inset_element(p,
                left = bg_coords[i, 1] - 0.04,
                right = bg_coords[i, 1] + 0.038,
                top = bg_coords[i, 2] + 0.04,
                bottom = bg_coords[i, 2] - 0.04,
                align_to = "panel"
            )
    }

    return(bg)
}


# Adapted from cowplot::theme_nothing
#' @importFrom ggplot2 '%+replace%'
# ggplot2 theme that strips all elements from a plot.
theme_blank <- function(title = ggplot2::element_blank(),
                        text = ggplot2::element_blank()) {
    ggplot2::theme_void() %+replace%
        ggplot2::theme(
            # Elements in this first block aren't used directly, but are inherited
            line = ggplot2::element_blank(),
            rect = ggplot2::element_rect(),
            text = text,
            aspect.ratio = 1,
            axis.line = ggplot2::element_blank(),
            axis.line.x = NULL,
            axis.line.y = NULL,
            axis.text = ggplot2::element_blank(),
            axis.text.x = NULL,
            axis.text.x.top = NULL,
            axis.text.y = NULL,
            axis.text.y.right = NULL,
            axis.ticks = ggplot2::element_blank(),
            axis.ticks.length = ggplot2::unit(0, "pt"),
            axis.title = ggplot2::element_blank(),
            axis.title.x = NULL,
            axis.title.x.top = NULL,
            axis.title.y = NULL,
            axis.title.y.right = NULL,
            legend.background = ggplot2::element_blank(),
            legend.spacing = NULL,
            legend.spacing.x = NULL,
            legend.spacing.y = NULL,
            legend.margin = ggplot2::margin(0, 0, 0, 0),
            legend.key = ggplot2::element_blank(),
            legend.key.size = NULL,
            legend.key.height = NULL,
            legend.key.width = NULL,
            legend.text = ggplot2::element_blank(),
            legend.text.align = NULL,
            legend.title = ggplot2::element_text(hjust = 0),
            legend.title.align = NULL,
            legend.position = "none",
            legend.direction = NULL,
            legend.justification = "center",
            legend.box = NULL,
            legend.box.margin = ggplot2::margin(0, 0, 0, 0),
            legend.box.background = ggplot2::element_blank(),
            legend.box.spacing = ggplot2::unit(0, "pt"),
            panel.grid = ggplot2::element_blank(),
            panel.grid.major = NULL,
            panel.grid.minor = NULL,
            panel.spacing = ggplot2::unit(0, "pt"),
            panel.spacing.x = NULL,
            panel.spacing.y = NULL,
            panel.ontop = FALSE,
            strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_blank(),
            strip.text.x = NULL,
            strip.text.y = NULL,
            strip.placement = "inside",
            strip.placement.x = NULL,
            strip.placement.y = NULL,
            strip.switch.pad.grid = ggplot2::unit(0., "cm"),
            strip.switch.pad.wrap = ggplot2::unit(0., "cm"),
            plot.background = ggplot2::element_blank(),
            plot.title = title,
            plot.subtitle = ggplot2::element_blank(),
            plot.caption = ggplot2::element_blank(),
            plot.tag = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(0, 0, 0, 0),
            panel.background = ggplot2::element_rect(
                fill = "#ffffffcc",
                colour = "#ffffffcc"
            ),
            panel.border = ggplot2::element_rect(
                colour = "black",
                fill = NA
            ),
            complete = TRUE
        )
}


# TODO: Improve color palette

#' ggplot2 scale for the MPI colors and extended swatches.
#' @param name The name of the scale. Either "mpi" or "mpimg".
#' @param ... Further arguments to ggplot2::discrete_scale.
scale_color_mpimg <- function(name = "mpimg", ...) {
    if (name == "mpimg") {
        ggplot2::discrete_scale(
            scale_name = "mpimg",
            aesthetics = "color",
            palette = mpimg_pal(),
            ...
        )
    } else if (name == "mpi") {
        ggplot2::discrete_scale(
            scale_name = "mpi",
            aesthetics = "color",
            palette = mpi_pal(),
            ...
        )
    }
}

#' MPIMG color palette.
#' @returns A function that can be used to generate colors.
mpimg_pal <- function() {
    mpi_colors <- c(
        "#006c66", # MPG-CD-Grün
        "#777777", # MPG-Dunkelgrau
        "#a7a7a8", # MPG-Grau
        "#c6d325", # MPG Hellgrün
        "#29485d", # MPG Dunkelblau
        "#00b1ea", # MPG Hellblau
        "#ef7c00" # MPG Orange
    )

    scales::manual_pal(values = mpi_colors)
}

#' MPI color palette
#' @returns A function that can be used to generate colors.
mpi_pal <- function() {
    # mpi colors extended.
    mpimg_colors <- c(
        "#006C66",
        "#29485D",
        "#009ACD",
        "#777777",
        "#C6D325",
        "#21AE2D",
        "#EF7C00",
        "#F5B742",
        "#57219D",
        "#950980"
    )

    scales::manual_pal(values = mpimg_colors)
}

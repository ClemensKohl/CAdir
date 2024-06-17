
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
cluster_apl_old <- function(caobj,
                        cadir,
                        cluster = NULL,
                        direction,
                        group,
                        show_cells = TRUE,
                        show_genes = FALSE,
                        show_lines = TRUE,
                        highlight_cluster = TRUE,
                        # colour_by_group = FALSE,
                        label_genes = FALSE,
                        point_size = 1.5,
                        size_factor = 2,
                        ntop = 15) {

    # TODO: Check if you can simplify the cluster/group assignments.
    stopifnot(methods::is(caobj, "cacomp"))
    stopifnot(methods::is(cadir, "cadir"))

    if (!is.null(cluster)) cluster <- as.character(cluster)

    all_cls <- as.character(
        sort(unique(c(cadir@cell_clusters, cadir@gene_clusters)))
    )

    all_cls <- c(levels(cadir@cell_clusters), levels(cadir@gene_clusters))

    if (!is.null(cluster) && (!cluster %in% all_cls)) cluster <- NULL
    if (is.null(cluster)) {
        # Kinda redundant. placeholder if I want to do deal with special case.
        cell_grp <- seq_len(length(cadir@cell_clusters))
        gene_grp <- seq_len(length(cadir@gene_clusters))
    } else {
        cell_grp <- which(cadir@cell_clusters == cluster)
        gene_grp <- which(cadir@gene_clusters == cluster)
    }

    # ensure that clusters and directions are coherent
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
    dapl_nms <- rownames(dapl)

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
            names(cadir@cell_clusters)[cadir@cell_clusters == dapl_nms[d]],
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
        # if (isTRUE(colour_by_group)) {
        #     sel <- group
        # } else {
            sel <- c()
            if (show_cells) {
                sel_cells <- match(
                    names(cadir@cell_clusters)[cell_grp],
                    df$sample
                )
                sel <- c(sel, sel_cells)
            }
            if (show_genes) {
                sel_genes <- match(
                    names(cadir@gene_clusters)[gene_grp],
                    df$sample
                )
                sel <- c(sel, sel_genes)
            }
        # }
        sel <- na.omit(sel)

        df$cluster <- "other"
        df$cluster[sel] <- "cluster"
        df$cluster <- factor(df$cluster, levels = c("other", "cluster"))
        ord <- order(df$cluster)
        df <- df[ord, ]
    } else if (show_cells || show_genes) {
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
            size_factor <- size_factor # I know, I know ...
        } else {
            size_factor <- 1
        }
        p <- p +
            ggplot2::geom_point(ggplot2::aes(shape = type, size = type, alpha = cluster)) +
            ggplot2::scale_size_manual(values = c(
                "cell" = point_size,
                "gene" = point_size * size_factor
            )) +
            ggplot2::scale_shape_manual(values = c(
                "cell" = 19,
                "gene" = 8
            )) +
            ggplot2::scale_alpha_manual(values = c(
                "cluster" = 1,
                "other" = 0.7,
                "cell_cluster" = 1,
                "gene_cluster" = 1,
                "cell_other" = 0.7,
                "gene_other" = 0.3
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

            if (label_genes) {
                to_highlight <- (df$type == "gene" &
                    df$cluster == "gene_cluster")
                dfh <- df[to_highlight, ]
                dfh <- head(dfh[order(dfh$x, decreasing = TRUE), ], ntop)
                p <- p + ggrepel::geom_label_repel(
                    data = dfh,
                    ggplot2::aes(
                        x = x,
                        y = y,
                        label = sample
                    ),
                    max.overlaps = Inf
                )
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
    #FIXME: Doesnt work with annotated data.
    size <- 1
    pls <- list()
    cls <- levels(cadir@cell_clusters)
    anno_dirs <- all(rownames(cadir@directions) %in% cls)

    for (i in seq_along(cls)) {
        for (j in seq_along(cls)) {
            sel <- which(cadir@cell_clusters == cls[i] |
                         cadir@cell_clusters == cls[j])

            if (anno_dirs) {
                sel_dir <- which(rownames(cadir@directions) %in% unique(c(cls[i], cls[j])))
                dir_idx <- which(rownames(cadir@directions) == cls[i])
            } else {
                sel_dir <- search_dict(cadir@dict, c(cls[i], cls[j]))
                dir_idx <- search_dict(cadir@dict, cls[i])
            }

            sub_cak <- methods::new("cadir",
                cell_clusters = cadir@cell_clusters[sel],
                directions = cadir@directions[sel_dir, , drop = FALSE]
            )
            # cat("\n\ni,j", i, ",", j)
            if (show_cells) sc <- i == j else sc <- FALSE
            if (show_genes) sg <- i == j else sg <- FALSE
            p <- cluster_apl(
                caobj = caobj,
                cadir = sub_cak,
                direction = cadir@directions[dir_idx, ],
                group = which(cadir@cell_clusters == cls[i]),
                cluster = cls[i],
                show_cells = sc,
                show_genes = sg,
                show_lines = i != j,
                highlight_cluster = highlight_cluster,
                # colour_by_group = TRUE
            ) +
                ggplot2::ggtitle(ifelse(i == j, paste0("Cluster: ", cls[i]), "")) +
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

    fig <- ggpubr::ggarrange(
        plotlist = pls,
        nrow = length(cls),
        ncol = length(cls)
    )
    return(fig)
}

#' Summarizes the cell clustering results in a single plot.
#' @param cadir A cadir object with valid cell clustering results.
#' @param caobj A cacomp object.
#' @param point_size Size of the points (cells).
#' @inheritParams cluster_apl
#' @returns A plot that summarizes the cell clustering results.
#' @export
plot_clusters <- function(cadir,
                          caobj,
                          point_size = 1,
                          size_factor = 1,
                          show_genes = FALSE,
                          label_genes = FALSE,
                          ntop = 5,
                          text_size = 16) {
    pls <- list()
    cls <- levels(cadir@cell_clusters)

    for (i in seq_along(cls)) {
        p <- cluster_apl(
            caobj = caobj,
            cadir = cadir,
            direction = cadir@directions[cls[i], ],
            cluster = as.character(cls[i]),
            group = which(cadir@cell_clusters == cls[i]),
            show_cells = TRUE,
            show_genes = show_genes,
            show_lines = FALSE,
            highlight_cluster = TRUE,
            # colour_by_group = FALSE,
            label_genes = label_genes,
            point_size = point_size,
            size_factor = size_factor,
            ntop = ntop
        ) +
            ggplot2::ggtitle(paste0("Cluster: ", cls[i])) +
            ggplot2::theme(
                legend.position = "none",
                axis.title.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                plot.title = ggplot2::element_text(size = text_size)
            )

        # p$layers[[1]]$aes_params$size <- point_size

        pls[[i]] <- p
    }

    fig <- ggpubr::ggarrange(
        plotlist = pls,
        nrow = ceiling(sqrt(length(cls))),
        ncol = ceiling(sqrt(length(cls)))
    )
    return(suppressWarnings(fig))
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
#' @export
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
                    org = "mm",
                    keep_end = TRUE) {

    # FIXME: Genes are basically impossible to tell from cells
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

            tmp_cadir@gene_clusters <- CAdir:::assign_genes(
                caobj = caobj,
                cadir = tmp_cadir,
                qcutoff = cadir@parameters$qcutoff
            )

            if (isTRUE(annotate_clusters)) {
                suppressWarnings({
                    tmp_cadir <- CAdir::annotate_biclustering(
                        obj = tmp_cadir,
                        universe = rownames(caobj@std_coords_rows),
                        org = org,
                        alpha = 0.05,
                        min_size = 10,
                        max_size = 500
                    )
                })
            }
            # else {
            #     #FIXME: After fixing direction naming, this should be redundant.
            #     tmp_cadir@cell_clusters <- factor(
            #         paste0("cluster_", f2c(tmp_cadir@cell_clusters)),
            #         levels = paste0("cluster_", levels(tmp_cadir@cell_clusters))
            #     )
            #     names(tmp_cadir@cell_clusters) <- names(tmp_ccs)
            #
            #     gnms <- names(cadir@gene_clusters)
            #     tmp_cadir@gene_clusters <- factor(
            #         paste0("cluster_", f2c(tmp_cadir@gene_clusters)),
            #         levels = paste0("cluster_", levels(tmp_cadir@gene_clusters))
            #     )
            #     names(tmp_cadir@gene_clusters) <- gnms
            #
            #     rownames(tmp_cadir@directions) <- gsub(
            #         pattern = "line",
            #         replacement = "cluster_",
            #         x = rownames(tmp_cadir@directions)
            #     )
            # }

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
            # colour_by_group = colour_by_group,
            show_lines = FALSE,
            point_size = 0.3
        )
        if (isTRUE(annotate_clusters)) {
            p <- p + ggplot2::ggtitle(cluster) +
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
#' @export
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

#' ggplot2 scale for the MPI colors and extended swatches.
#' @param name The name of the scale. Either "mpi" or "mpimg".
#' @param ... Further arguments to ggplot2::discrete_scale.
#' @export
scale_fill_mpimg <- function(name = "mpimg", ...) {
    if (name == "mpimg") {
        ggplot2::discrete_scale(
            scale_name = "mpimg",
            aesthetics = "fill",
            palette = mpimg_pal(),
            ...
        )
    } else if (name == "mpi") {
        ggplot2::discrete_scale(
            scale_name = "mpi",
            aesthetics = "fill",
            palette = mpi_pal(),
            ...
        )
    }
}

#' MPIMG color palette.
#' @returns A function that can be used to generate colors.
#' @export
mpimg_pal <- function() {
    mpi_colors <- c(
        "#006c66", # MPG-CD-Grün
        "#c6d325", # MPG Hellgrün
        "#ef7c00", # MPG Orange
        "#29485d", # MPG Dunkelblau
        "#00b1ea", # MPG Hellblau
        "#777777", # MPG-Dunkelgrau
        "#a7a7a8" # MPG-Grau
    )

    scales::manual_pal(values = mpi_colors)
}

#' MPI color palette
#' @returns A function that can be used to generate colors.
#' @export
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

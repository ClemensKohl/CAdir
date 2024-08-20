#' Plot for the clustering results wich shows the
#' relationship between the clusters. In the diagonal
#' an APL plot for the respective cluster is shown.
#' @param cadir A cadir object with valid cell clustering results.
#' @param caobj A cacomp object.
#' @inheritParams cluster_apl
#' @inheritParams plot_clusters
#' @param ... Further arguments forwarded to cluster_apl().
#' @returns A plot that summarizes the cell clustering results and
#' how the clusters relate to each other.
#' @export
plot_results <- function(cadir,
                         caobj,
                         highlight_cluster = TRUE,
                         show_cells = TRUE,
                         show_genes = FALSE,
                         title_prefix = "Cluster: ",
                         gsub_title = "",
                         return_list = FALSE,
                         ...
                         ) {
    base::stopifnot(
        "Set either `show_cells` or `show_genes` to TRUE." =
            isTRUE(show_cells) || isTRUE(show_genes)
    )
    size <- 1
    pls <- list()
    cls <- levels(cadir@cell_clusters)
    anno_dirs <- all(rownames(cadir@directions) %in% cls)

    for (i in seq_along(cls)) {
        for (j in seq_along(cls)) {

            if (!is.null(gsub_title)) {
                cls_title <- gsub(
                                  pattern = gsub_title,
                                  replacement = " ",
                                  x = cls[i]
                )
            } else {
                cls_title <- cls[i]
            }

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
                ...
            ) +
                ggplot2::ggtitle(ifelse(i == j, paste0(title_prefix, cls_title), "")) +
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

    if (isTRUE(return_list)) {
        fig <- pls
    } else {
        fig <- ggpubr::ggarrange(
            plotlist = pls,
            nrow = length(cls),
            ncol = length(cls)
    )

    }

    return(fig)
}

#' Summarizes the cell clustering results in a single plot.
#' @param cadir A cadir object with valid cell clustering results.
#' @param caobj A cacomp object.
#' @param point_size Size of the points (cells).
#' @param text_size Size of the text in the plot.
#' @param title_prefix Prefix to print before cluster name.
#' @param ggncol Number of columns to arrange plots in.
#' @param ggnrow Number of rows to arrange plots in.
#' @param axis Whether to show axis markings or not.
#' @param gsub_title Character string to substitute with " " in plot titles.
#' @param legend_pos Where to position the legend.
#' Normal ggplot positions allowed.
#' @param return_list If TRUE, instead of a final, arranged panel the plots
#' are returned as a list.
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
                          text_size = 16,
                          title_prefix = "Cluster: ",
                          ggncol = NULL,
                          ggnrow = NULL,
                          axis = FALSE,
                          gsub_title = NULL,
                          legend_pos = "none",
                          return_list = FALSE) {
    pls <- list()
    cls <- levels(cadir@cell_clusters)

    if (isFALSE(axis)) {
        plot_theme <- ggplot2::theme(
            legend.position = legend_pos,
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = text_size)
        )
    } else {
        plot_theme <- ggplot2::theme(
            legend.position = legend_pos,
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = text_size)
        )
    }
    for (i in seq_along(cls)) {
        if (!is.null(gsub_title)) {
            cls_title <- gsub(
                              pattern = gsub_title,
                              replacement = " ",
                              x = cls[i]
            )
        } else {
            cls_title <- cls[i]
        }
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
            ggplot2::ggtitle(paste0(title_prefix, cls_title)) +
            plot_theme

        pls[[i]] <- p
    }

    if (isTRUE(return_list)) {
        fig <- pls
    } else {
        fig <- ggpubr::ggarrange(
            plotlist = pls,
            nrow = ifelse(test = is.null(ggnrow),
                        yes = ceiling(sqrt(length(cls))),
                        no = ggnrow),
            ncol = ifelse(test = is.null(ggncol),
                        yes = ceiling(sqrt(length(cls))),
                        no = ggncol)
        )
    }

    return(suppressWarnings(fig))
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

#' @importFrom ggplot2 '%+replace%'
# ggplot2 theme that strips all elements from a plot, but leaves axis ticks and elements.
theme_axis_only <- function(title = ggplot2::element_blank(),
                        text = ggplot2::element_blank()) {
    ggplot2::theme_classic() %+replace%
        ggplot2::theme(
            # Elements in this first block aren't used directly, but are inherited
            # line = ggplot2::element_blank(),
            rect = ggplot2::element_rect(),
            text = text,
            aspect.ratio = 1,
            axis.line = ggplot2::element_blank(),
            axis.line.x = NULL,
            axis.line.y = NULL,
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
    } else if (name == "mpi_extend") {
        ggplot2::discrete_scale(
            scale_name = "mpi_extend",
            aesthetics = "color",
            palette = mpi_extend_pal(),
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
    } else if (name == "mpi_extend") {
        ggplot2::discrete_scale(
            scale_name = "mpi_extend",
            aesthetics = "fill",
            palette = mpi_extend_pal(),
            ...
        )
    }
}


scale_fill_gradient_mpimg <- function(name = "orange", ...) {
    if (name == "orange") {
        ggplot2::scale_fill_gradientn(
            colours = c("#29485d", "#006c66", "#c6d325", "#ef7c00"),
            ...
        )
    } else if (name == "green") {
        ggplot2::scale_fill_gradientn(
            colours = c("#29485d", "#006c66", "#c6d325"),
            ...
        )
    } else {
        stop("Pick a correct name!")
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
mpi_extend_pal <- function() {
    # mpi colors extended.
    mpi_extend_colors <- c(
        "#006c66", # MPG-CD-Grün
        "#c6d325", # MPG Hellgrün
        "#ef7c00", # MPG Orange
        "#29485d", # MPG Dunkelblau
        "#00b1ea", # MPG Hellblau
        "#777777", # MPG-Dunkelgrau
        "#d44a3d", # Warm Red
        "#8c5fa8", # Soft Purple
        "#a7a7a8", # MPG-Grau
        "#f4c542", # Golden Yellow
        "#d95276", # Rich Pink
        "#a25b43", # Warm Brown
        "#00CED1" # Soft Cyan
    )

    scales::manual_pal(values = mpi_extend_colors)
}

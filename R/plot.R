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
    group = NULL
) {
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
#' @param cluster_id The cluster for which to plot the APL.
#' @param cluster_id The name of the cluster for the plot title.
#' @param show_points If TRUE, points (cells) are plotted.
#' @param show_lines If TRUE, the directions in cadir are plotted.
#' @param plot_group If TRUE, highlights the points in `group`.
#' @param point_size Size of the points (cells).
#' @returns
#' An APL plot (ggplot2 object).
cluster_apl <- function(
        caobj,
        cadir,
        direction,
        group,
        cluster_id = "NA",
        show_points = TRUE,
        show_lines = TRUE,
        plot_group = FALSE,
        point_size = 1.5) {

    stopifnot(methods::is(caobj, "cacomp"))
    stopifnot(methods::is(cadir, "cadir"))

    # ensure that clusters and directions are coherent
    cadir <- rename_clusters(cadir)
    if (nrow(cadir@directions) > 1) {

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

    capl <- model(caobj@prin_coords_cols)
    dapl <- model(cadir@directions)

    # If the line points into the opposite direction of points
    # we flip the line.
    for (d in seq_len(nrow(dapl))) {
        sel <- match(
            names(cadir@cell_clusters)[f2n(cadir@cell_clusters) == d],
            rownames(caobj@prin_coords_cols)
        )

        if (length(sel) == 0) next

        if (length(sel) > 1) {
            grp_mean <- colMeans(capl[sel, ])
        } else {
            grp_mean <- capl[sel, ]
        }

        if (sign(grp_mean[1]) != sign(dapl[d, 1])) {
            dapl[d, ] <- c(-1, 1) * dapl[d, ]
        }
    }

    df <- as.data.frame(capl)
    df$sample <- rownames(df)

    if (isTRUE(plot_group)) {
        sel <- match(rownames(caobj@prin_coords_cols)[group], df$sample)
        df$cluster <- "other"
        df$cluster[sel] <- "cluster"
    } else {
        df$cluster <- 0
        sel <- match(names(cadir@cell_clusters), df$sample)
        df$cluster[sel] <- cadir@cell_clusters
    }

    df$cluster <- as.factor(df$cluster)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = cluster))

    if (isTRUE(show_points)) {

        p <- p +
            ggplot2::geom_point(size = point_size)

        if (isTRUE(plot_group)) {
            p <- p + ggplot2::scale_color_manual(values = c("cluster" = "#c6d325",
                                                            "other" = "#006c66"))
        }
    }

    if (isTRUE(show_lines)) {

        for (d in seq_len(nrow(dapl))) {

            if (isTRUE(all.equal(dapl[d, ], c(1, 0), tolerance = 1e-4, check.attributes = FALSE)) ||
                isTRUE(all.equal(dapl[d, ], c(-1, 0), tolerance = 1e-4, check.attributes = FALSE))) {
                lcolor <- "black"
                ltype <- "solid"
            } else {
                lcolor <- "red"
                ltype <- "dashed"

                if(isTRUE(plot_group)) {
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
            cluster_id,
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
#' @returns A plot that summarizes the cell clustering results and
#' how the clusters relate to each other.
plot_results <- function(cadir, caobj) {

  size <- 1
  pls <- list()
  cls <- sort(unique(cadir@cell_clusters))

  for (i in seq_along(cls)) {

    for (j in seq_along(cls)) {

      sel <- which(cadir@cell_clusters == cls[i] | cadir@cell_clusters == cls[j])
      sel_dir <- unique(c(f2n(cls[i]), f2n(cls[j])))

      sub_cak <- methods::new("cadir",
                     cell_clusters = cadir@cell_clusters[sel],
                     directions = cadir@directions[sel_dir, , drop = FALSE])

      p <- cluster_apl(caobj = caobj,
                       cadir = sub_cak,
                       direction = cadir@directions[f2n(cls[i]), ],
                       group = which(cadir@cell_clusters == cls[i]),
                       cluster_id = "cell_clusters",
                       show_points = i == j,
                       show_lines = i != j,
                       plot_group = TRUE) +
                           ggplot2::ggtitle("") +
                           ggplot2::theme(legend.position = "none",
                                 axis.title.x = ggplot2::element_blank(),
                                 axis.text.x = ggplot2::element_blank(),
                                 axis.ticks.x = ggplot2::element_blank(),
                                 axis.title.y = ggplot2::element_blank(),
                                 axis.text.y = ggplot2::element_blank(),
                                 axis.ticks.y = ggplot2::element_blank())


                           if (i == j) {
                             p$layers[[1]]$aes_params$size <- size
                           } else {
                             p$layers[[2]]$aes_params$size <- size
                           }
                           pls[[paste0(i, "_", j)]] <- p
    }

  }

  fig <- ggpubr::ggarrange(plotlist = pls, nrow = length(cls), ncol = length(cls))
  return(fig)
}

#' Summarizes the cell clustering results in a single plot.
#' @param cadir A cadir object with valid cell clustering results.
#' @param caobj A cacomp object.
#' @param point_size Size of the points (cells).
#' @returns A plot that summarizes the cell clustering results.
plot_clusters <- function(cadir, caobj, point_size = 1) {
    pls <- list()
    cls <- sort(unique(cadir@cell_clusters))

    for (i in seq_along(cls)) {

        p <- cluster_apl(caobj = caobj,
                         cadir = cadir,
                         direction = cadir@directions[f2n(cls[i]), ],
                         group = which(cadir@cell_clusters == cls[i]),
                         show_points = TRUE,
                         show_lines = FALSE,
                         plot_group = TRUE) +
                      ggplot2::ggtitle(paste0("cluster_", i)) +
                      ggplot2::theme(legend.position = "none",
                                     axis.title.x = ggplot2::element_blank(),
                                     axis.text.x = ggplot2::element_blank(),
                                     axis.ticks.x = ggplot2::element_blank(),
                                     axis.title.y = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_blank(),
                                     axis.ticks.y = ggplot2::element_blank())

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


# TODO: Move to a package that is on CRAN (so not ggsankey)

#' Creates a sankey plot of the clustering results with all splits and merges.
#' @param cadir A cadir object with valid clustering results.
#' @param rm_redund If TRUE, only shows an iteration if something changes.
#' @returns A sankey plot of the clustering results.
plot_flow <- function(cadir, rm_redund = TRUE) {
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
        sub_cls <- sub_cls[, -del_cl]
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

plot_networkd3 <- function(cadir, rm_redund = rm_redund) {

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
        sub_cls <- sub_cls[, -del_cl]
    }
    nodes <- data.frame("name" = unique(as.vector(as.matrix(sub_cls))))
    links <- data.frame()
    for (n in seq_len(ncol(sub_cls))) {
        from  <- unique(sub_cls[, n])

        for (l in seq_len(length(from))) {
            idx <- which(sub_cls[, n] == from[l])
            to <- unique(sub_cls[idx, n + 1])

            for (t in seq_len(length(to))) {
                val <- sum(sub_cls[idx, n + 1] == to[t])

                tmpdf <- data.frame(
                    "source" = from[l],
                    "target" = to[t],
                    "value" = val
                )

                links <- rbind(links, tmpdf)
            }
        }

    }



    sankeyNetwork(
        Links = links,
        Nodes = nodes,
        Source = "source",
        Target = "target",
        Value = "value",
        NodeID = "name",
        fontSize = 12,
        nodeWidth = 30
    )



}

#' Plots the graph of the clustering splits and merges
#' and overlays APL plots over the nodes.
#' @inheritParams plot_sm_graph
#' @param caobj A cacomp object.
#' @returns
#' A ggplot object showing the split-merge graph and APL plots for each cluster.
sm_plot <- function(cadir, caobj, rm_redund = TRUE) {

    graph <- build_graph(cadir = cadir, rm_redund = rm_redund)

    lgraph <- ggraph::create_layout(graph, layout="tree")

    ggraph::set_graph_style(plot_margin = ggplot2::margin(0, 0, 0, 0))
    bg <- ggraph::ggraph(lgraph) +
        ggraph::geom_edge_link() +
        ggraph::geom_node_point(alpha = 1)

    bg_coords <- get_x_y_values(bg)

    cls <- cadir@log$clusters

    nodes <- names(igraph::V(graph))

    for (i in seq_len(nrow(lgraph))) {

        node_nm <- nodes[i]
        # TODO: replace with base R code and remove dependency.
        name_elems <- stringr::str_split_1(node_nm, "-")

        if (name_elems[1] == "root") next

        iter_nm  <- name_elems[1]
        cluster <- as.numeric(name_elems[2])

        grp_idx <- which(cls[, iter_nm] == cluster)

        # TODO: We calculate the directions new. Shouldn't we save them somehow?
        dir <- total_least_squares(caobj@prin_coords_cols[grp_idx, ])

        p <- cluster_apl(caobj = caobj,
                         cadir = cadir,
                         direction = as.numeric(dir),
                         group = grp_idx,
                         cluster_id = as.character("cluster"),
                         show_lines = FALSE,
                         point_size = 0.3) +
                      scale_color_mpimg(name = "mpimg") +
                      theme_blank()

                      bg <-  bg +
                          patchwork::inset_element(p,
                                                   left = bg_coords[i, 1] - 0.04,
                                                   right = bg_coords[i, 1] + 0.038,
                                                   top = bg_coords[i, 2] + 0.04,
                                                   bottom = bg_coords[i, 2] - 0.04,
                                                   align_to = "panel")

    }

    return(bg)

}


# Adapted from cowplot::theme_nothing
#' @importFrom ggplot2 '%+replace%'
# ggplot2 theme that strips all elements from a plot.
theme_blank <- function() {
    ggplot2::theme_void() %+replace%
        ggplot2::theme(
              # Elements in this first block aren't used directly, but are inherited
              line = ggplot2::element_blank(),
              rect = ggplot2::element_rect(),
              text = ggplot2::element_blank(),
              aspect.ratio = 1,
              axis.line =          ggplot2::element_blank(),
              axis.line.x =        NULL,
              axis.line.y =        NULL,
              axis.text =          ggplot2::element_blank(),
              axis.text.x =        NULL,
              axis.text.x.top =    NULL,
              axis.text.y =        NULL,
              axis.text.y.right =  NULL,
              axis.ticks =         ggplot2::element_blank(),
              axis.ticks.length =  ggplot2::unit(0, "pt"),
              axis.title =         ggplot2::element_blank(),
              axis.title.x =       NULL,
              axis.title.x.top =   NULL,
              axis.title.y =       NULL,
              axis.title.y.right = NULL,

              legend.background =  ggplot2::element_blank(),
              legend.spacing =     NULL,
              legend.spacing.x =   NULL,
              legend.spacing.y =   NULL,
              legend.margin =      ggplot2::margin(0, 0, 0, 0),
              legend.key =         ggplot2::element_blank(),
              legend.key.size =   NULL,
              legend.key.height =  NULL,
              legend.key.width =   NULL,
              legend.text =        ggplot2::element_blank(),
              legend.text.align =  NULL,
              legend.title =       ggplot2::element_text(hjust = 0),
              legend.title.align = NULL,
              legend.position =    "none",
              legend.direction =   NULL,
              legend.justification = "center",
              legend.box =         NULL,
              legend.box.margin =  ggplot2::margin(0, 0, 0, 0),
              legend.box.background = ggplot2::element_blank(),
              legend.box.spacing = ggplot2::unit(0, "pt"),

              panel.grid =         ggplot2::element_blank(),
              panel.grid.major =   NULL,
              panel.grid.minor =   NULL,
              panel.spacing =      ggplot2::unit(0, "pt"),
              panel.spacing.x =    NULL,
              panel.spacing.y =    NULL,
              panel.ontop    =     FALSE,

              strip.background =   ggplot2::element_blank(),
              strip.text =         ggplot2::element_blank(),
              strip.text.x =       NULL,
              strip.text.y =       NULL,
              strip.placement =    "inside",
              strip.placement.x =  NULL,
              strip.placement.y =  NULL,
              strip.switch.pad.grid = ggplot2::unit(0., "cm"),
              strip.switch.pad.wrap = ggplot2::unit(0., "cm"),

              plot.background =    ggplot2::element_blank(),
              plot.title =         ggplot2::element_blank(),
              plot.subtitle =      ggplot2::element_blank(),
              plot.caption =       ggplot2::element_blank(),
              plot.tag           = ggplot2::element_blank(),
              plot.margin =        ggplot2::margin(0, 0, 0, 0),

              panel.background = ggplot2::element_rect(fill = "#ffffffcc",
                                                       colour = "#ffffffcc"),
              panel.border = ggplot2::element_rect(colour = "black",
                                                   fill = NA),
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

# TODO: add documentation
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

# TODO: add documentation
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

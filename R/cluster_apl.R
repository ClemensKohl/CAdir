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
#' @param point_size Size of the points (cells).
#' @param size_factor Factor by which the genes are
#' larger/smaller than `point_size`.
#' @param label_genes If TRUE, adds text labels for the
#' `ntop` genes per cluster.
#' @param ntop Number of genes to label if `label_genes = TRUE`.
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
                        label_genes = FALSE,
                        point_size = 1.5,
                        size_factor = 2,
                        ntop = 15) {
    stopifnot(methods::is(caobj, "cacomp"))
    stopifnot(methods::is(cadir, "cadir"))

    if (is.null(cluster) && isTRUE(highlight_cluster)) {
        warning(paste0("Turning `highlight_cluster off,",
                       "because no `cluster` was specified"))

        highlight_cluster <- FALSE
    }

    all_cls <- unique(c(
        levels(cadir@cell_clusters),
        levels(cadir@gene_clusters)
    ))

    if (!is.null(cluster)) {
        cluster <- as.character(cluster)
        if (!cluster %in% all_cls) cluster <- NULL
    }


    # Calculate angle if only two directions, NA otherwise
    ang <- .get_plot_angle(directions = cadir@directions)

    model <- apl_model(
        caobj = caobj,
        direction = direction,
        group = group
    )

    dapl <- model(cadir@directions)
    dapl <- .toggle_dir(
        apl_dirs = dapl,
        caobj = caobj,
        cadir = cadir,
        model = model
    )

    bool_sum <- show_cells + show_genes
    coords <- list("prin_coords_cols", "std_coords_cols")

    c_coords <- if (isTRUE(show_cells)) {
        slot(caobj, name = coords[[bool_sum]])
    } else {
        NULL
    }

    df <- .construct_df(
        cells = c_coords,
        genes = if (isTRUE(show_genes)) caobj@prin_coords_rows else NULL,
        model = model
    )

    df <- .add_cluster_info(
        plot_df = df,
        cadir = cadir,
        cluster = cluster,
        highlight_cluster = highlight_cluster
    )

    ############
    ### Plot ###
    ############

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = cluster))

    p <- .cluster_apl_points(
        ggplt = p,
        plot_df = df,
        point_size = point_size,
        size_factor = size_factor,
        highlight_cluster = highlight_cluster,
        label_genes = label_genes,
        ntop = ntop
    )

    if (isTRUE(show_lines)) {
        p <- .add_lines(
            ggplt = p,
            apl_dir = dapl,
            highlight_cluster = highlight_cluster
        )
    }

    p <- p +
        ggplot2::ggtitle(paste0(
            "Cluster: ",
            as.character(cluster),
            ", CA-angle: ",
            round(ang, 2)
        )) +
        ggplot2::theme_bw()

    return(p)
}


# TODO: Add documentation
.construct_df <- function(cells = NULL, genes = NULL, model = NULL) {
    incl_cells <- !is.null(cells)
    incl_genes <- !is.null(genes)

    if (incl_cells && !incl_genes) {
        capl <- model(cells)
        df <- as.data.frame(capl)
        df$sample <- rownames(df)
        df$type <- "cell"
    } else if (incl_cells && incl_genes) {
        capl <- model(cells)
        gapl <- model(genes)
        df <- as.data.frame(capl)
        df$sample <- rownames(df)
        df$type <- "cell"

        dfg <- as.data.frame(gapl)
        dfg$sample <- rownames(dfg)
        dfg$type <- "gene"

        df <- rbind(df, dfg)
    } else if (!incl_cells && incl_genes) {
        gapl <- model(genes)
        df <- as.data.frame(gapl)
        df$sample <- rownames(df)
        df$type <- "gene"
    } else {
        df <- data.frame(
            "sample" = c(),
            "type" = c()
        )
    }

    return(df)
}

# TODO: Add documentation.
.get_plot_angle <- function(directions) {
    if (nrow(directions) == 2) {
        ang <- min(
            rad2deg(get_angle(directions[1, ], directions[2, ])),
            rad2deg(get_angle(-directions[1, ], directions[2, ]))
        )
    } else {
        ang <- NA_real_
    }

    return(ang)
}

# TODO: Add documentation.
.toggle_dir <- function(apl_dirs, caobj, cadir, model) {
    # If the line points into the opposite direction of points
    # we flip the line.
    dapl_nms <- rownames(apl_dirs)
    for (d in seq_len(nrow(apl_dirs))) {
        sel <- match(
            names(cadir@cell_clusters)[cadir@cell_clusters == dapl_nms[d]],
            rownames(caobj@prin_coords_cols)
        )

        if (length(sel) == 0) next

        cell_coords <- model(caobj@prin_coords_cols)[sel, ]

        if (length(sel) > 1) {
            grp_mean <- colMeans(cell_coords)
        } else {
            grp_mean <- cell_coords
        }

        if (sign(grp_mean[1]) != sign(apl_dirs[d, 1])) {
            apl_dirs[d, ] <- c(-1, 1) * apl_dirs[d, ]
        }
    }

    return(apl_dirs)
}

# TODO: Add documentation.
.add_cluster_info <- function(
    plot_df,
    cadir,
    cluster = NULL,
    highlight_cluster = FALSE) {
    show_genes <- any("gene" %in% plot_df$type)
    show_cells <- any("cell" %in% plot_df$type)

    if (!any(c(show_genes, show_cells))) {
        return(plot_df)
    }

    highlight_cluster <- !is.null(cluster) && isTRUE(highlight_cluster)
    if (isTRUE(highlight_cluster)) {
        cell_cl <- which(cadir@cell_clusters == cluster)
        gene_cl <- which(cadir@gene_clusters == cluster)
    } else {
        cell_cl <- seq_len(length(cadir@cell_clusters))
        gene_cl <- seq_len(length(cadir@gene_clusters))
    }

    sel_cells <- match(names(cadir@cell_clusters)[cell_cl], plot_df$sample)
    sel_cells <- na.omit(sel_cells)
    sel_genes <- match(names(cadir@gene_clusters)[gene_cl], plot_df$sample)
    sel_genes <- na.omit(sel_genes)


    if (isTRUE(highlight_cluster)) {
        sel <- c() # empty in case we only want to plot lines.
        sel <- c(sel, sel_cells, sel_genes)
        sel <- na.omit(sel)

        plot_df$cluster <- "other"
        plot_df$cluster[sel] <- "cluster"
        plot_df$cluster <- factor(plot_df$cluster,
            levels = c("other", "cluster")
        )

        ord <- order(plot_df$cluster)
        plot_df <- plot_df[ord, ]

        # if conditions, conc. type and cluster.
        if (show_cells && show_genes) {
            plot_df$cluster <- paste0(plot_df$type, "_", plot_df$cluster)
        }
    } else {
        no_cluster <- "N/A"
        plot_df$cluster <- no_cluster

        plot_df$cluster[sel_cells] <- f2c(cadir@cell_clusters)
        plot_df$cluster[sel_genes] <- f2c(cadir@gene_clusters)

        plot_df$cluster <- factor(plot_df$cluster,
            levels = sort(unique(c(
                no_cluster,
                f2c(cadir@cell_clusters),
                f2c(cadir@gene_clusters)
            )))
        )
    }

    return(plot_df)
}

# TODO: Add documentation.
.cluster_apl_points <- function(
    ggplt,
    plot_df,
    point_size = 1.5,
    size_factor = 1,
    highlight_cluster = FALSE,
    label_genes = FALSE,
    ntop = 15) {
    show_genes <- any("gene" %in% plot_df$type)
    show_cells <- any("cell" %in% plot_df$type)


    if (isTRUE(show_cells) || isTRUE(show_genes)) {
        ggplt <- ggplt +
            ggplot2::geom_point(ggplot2::aes(
                shape = type,
                size = type,
            ), alpha = 0.7)

        ggplt <- .add_colors(
            ggplt = ggplt,
            point_size = point_size,
            size_factor = size_factor
        )

        if (isTRUE(highlight_cluster)) {
            ggplt <- .highlight_cluster(
                ggplt = ggplt,
                both = (show_cells && show_genes)
            )

            if (isTRUE(label_genes)) {
                ggplt <- .label_genes(
                    ggplt = ggplt,
                    plot_df = plot_df,
                    ntop = ntop
                )
            }
        }
    }

    return(ggplt)
}

# TODO: Add documentation.
.add_colors <- function(ggplt, point_size = 1.5, size_factor = 2) {
    ggplt <- ggplt +
        ggplot2::scale_size_manual(values = c(
            "cell" = point_size,
            "gene" = point_size * size_factor
        )) +
        ggplot2::scale_shape_manual(values = c(
            "cell" = 19,
            "gene" = 8
        ))
    return(ggplt)
}

# TODO: Add documentation.
.highlight_cluster <- function(ggplt, both = FALSE) {
    if (isTRUE(both)) {
        ggplt <- ggplt + ggplot2::scale_color_manual(values = c(
            "cell_cluster" = "#c6d325",
            "gene_cluster" = "#ef7c00",
            "cell_other" = "#006c66",
            "gene_other" = "#777777"
        )) # +
        # ggplot2::aes(alpha = cluster) +
        # ggplot2::scale_alpha_manual(values = c(
        #     "cell_cluster" = 1,
        #     "gene_cluster" = 1,
        #     "cell_other" = 0.7,
        #     "gene_other" = 0.3
        # ))
    } else {
        ggplt <- ggplt +
          ggplot2::scale_color_manual(values = c(
            "cluster" = "#c6d325",
            "other" = "#006c66"
            ))
            # +
            # ggplot2::aes(alpha = cluster) +
            # ggplot2::scale_alpha_manual(values = c(
            #     "cluster" = 1,
            #     "other" = 0.7,
            # ))
    }

    return(ggplt)
}


# TODO: Add documentation.
.label_genes <- function(ggplt, plot_df, ntop = 15) {
    to_highlight <- (plot_df$cluster == "gene_cluster")

    dfh <- plot_df[to_highlight, ]
    dfh <- head(dfh[order(dfh$x, decreasing = TRUE), ], ntop)
    ggplt <- ggplt + ggrepel::geom_label_repel(
        data = dfh,
        ggplot2::aes(
            x = x,
            y = y,
            label = sample
        ),
        max.overlaps = Inf
    )

    return(ggplt)
}

# TODO: Add documentation.
.add_lines <- function(ggplt, apl_dir, highlight_cluster = FALSE) {
    for (d in seq_len(nrow(apl_dir))) {
        is_x <- is_xaxis(apl_dir[d, ])

        if (is_x) {
            lcolor <- "black"
            ltype <- "solid"
        } else {
            lcolor <- "red"
            ltype <- "dashed"
            if (isTRUE(highlight_cluster)) {
                lcolor <- "#006c66"
            }
        }

        ggplt <- ggplt + ggplot2::geom_abline(
            intercept = 0,
            slope = slope(lines = apl_dir[d, ], dims = 1:2),
            color = lcolor,
            linetype = ltype,
            size = 1
        ) +
            ggplot2::geom_point(
                data = data.frame(x = 0, y = 0),
                ggplot2::aes(x, y),
                color = lcolor
            )
    }

    return(ggplt)
}

# TODO: Add documentation.
is_xaxis <- function(apl_dir) {
    pos_xaxis <- all.equal(
        apl_dir,
        c(1, 0),
        tolerance = 1e-4,
        check.attributes = FALSE
    )

    neg_xaxis <- all.equal(
        apl_dir,
        c(-1, 0),
        tolerance = 1e-4,
        check.attributes = FALSE
    )

    return(isTRUE(pos_xaxis) || isTRUE(neg_xaxis))
}

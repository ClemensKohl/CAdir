#' Plot a cluster with the respective direction/line in an APL.
#' @param cadir A cadir object for which to compute the APL
#' @param caobj A cacomp object.
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
#' @param interactive Uses plotly to generate an interactive version of the
#' plot.
#' @returns
#' An APL plot (ggplot2 object).
#' @export
cluster_apl <- function(caobj,
                        cadir,
                        cluster,
                        direction = cadir@directions[cluster, ],
                        group = which(cadir@cell_clusters == cluster),
                        show_cells = TRUE,
                        show_genes = FALSE,
                        show_lines = FALSE,
                        highlight_cluster = TRUE,
                        label_genes = FALSE,
                        point_size = 1.5,
                        size_factor = 2,
                        ntop = 10,
                        interactive = FALSE) {
    stopifnot(methods::is(caobj, "cacomp"))
    stopifnot(methods::is(cadir, "cadir"))


    if (is.null(cluster) && isTRUE(highlight_cluster)) {
        rlang::warn(paste0("Turning `highlight_cluster off,",
                           "because no `cluster` was specified"))

        highlight_cluster <- FALSE
    }

    # Check if cluster is in actual cluster names.
    all_cls <- unique(c(
        levels(cadir@cell_clusters),
        levels(cadir@gene_clusters)
    ))

    if (!is.null(cluster)) {
        cluster <- as.character(cluster)
        if (!cluster %in% all_cls) {
            cluster <- NULL
            rlang::warn(paste0(
                "`cluster` does not correspond to a cluster in cadir.",
                " Setting to NULL."
            ))
        }
    }

    # Check validity of direction.
    ndim <- ncol(cadir@directions)
    dir_is_vec <- is.numeric(direction)
    has_len <- length(direction) == ndim
    has_grp <- !(length(group) == 0 || is.null(group))

    # Check if group is available.
    if (isTRUE(has_grp) && is.null(cluster)) {
        if (isFALSE(dir_is_vec) || isFALSE(has_len)){
            direction <- colMeans(caobj@prin_coords_cols[group, ])
            dir_is_vec <- is.numeric(direction)
            has_len <- length(direction) == ndim
        }
    } else if (isFALSE(has_grp)){
        rlang::abort("`group` has length 0.")
    }

    if (!dir_is_vec || !has_len) {
        rlang::abort("A valid direction is needed!")
    }


    # get gene ranks if we want to label them.
    if (isTRUE(label_genes) && length(cadir@gene_ranks) == 0) {
        cadir <- rank_genes(cadir = cadir, caobj = caobj)
    }

    # Calculate angle if only two directions, NA otherwise
    ang <- .get_plot_angle(directions = cadir@directions)

    bool_sum <- show_cells + show_genes
    coords <- list("prin_coords_cols", "std_coords_cols")

    coords_type <- if (bool_sum != 0) {
        coords[[bool_sum]]
    } else {
        "prin_coords_cols"
    }

    if (coords_type == "std_coords_cols" && !is.null(coords_type)) {
        # Convert direction to standard coordinates
        direction <- direction / caobj@D
    }

    model <- apl_model(
        coords = methods::slot(caobj, name = coords_type),
        direction = direction,
        group = group
    )


    c_coords <- if (isTRUE(show_cells)) {
        methods::slot(caobj, name = coords_type)
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
        highlight_cluster = highlight_cluster,
        label_genes = label_genes
    )

    if (nrow(df) > 0) {
        df$info <- paste0(
            "Name: ", df$sample, "\n",
            "Cluster: ", df$cluster, "\n",
            "Type: ", df$type
        )
    }

    ############
    ### Plot ###
    ############

    p <- ggplot2::ggplot(df, ggplot2::aes(
        x = x,
        y = y,
        color = cluster,
        text = info
    ))

    p <- .cluster_apl_points(
        ggplt = p,
        plot_df = df,
        point_size = point_size,
        size_factor = size_factor,
        highlight_cluster = highlight_cluster,
        label_genes = label_genes,
        ntop = ntop
    )

    if (!isFALSE(show_lines)) {
        show_both <- isTRUE(show_cells && show_genes)

        dapl <- cadir@directions
        if (coords_type == "std_coords_cols") {
            # Convert directions to standard coordinates
            dapl <- sweep(
                dapl,
                2,
                caobj@D,
                "/"
            )
        }

        dapl <- model(dapl)

        dapl <- .toggle_dir(
            apl_dirs = dapl,
            caobj = caobj,
            cadir = cadir,
            model = model,
            coords_type = coords_type
        )

        p <- .add_lines(
            ggplt = p,
            apl_dir = dapl,
            highlight_cluster = highlight_cluster,
            show_lines = show_lines,
            show_both = show_both
        )
    }

    p <- p +
        ggplot2::ggtitle(paste0(
            "Cluster: ",
           as.character(cluster)
            # ", CA-angle: ",
            # round(ang, 2)
        )) +
        ggplot2::theme_bw()

    if (isTRUE(interactive)) {
        suppressWarnings({
            p <- plotly::ggplotly(p = p, tooltip = "text")
        })
    }
    return(p)
}


#' Contructs the basic plotting data frame for cluster_apl
#' @param cells Cell coordinates. If NULL cells will be omitted.
#' @param genes Gene coordinates. If NULL genes will be omitted.
#' @param model The APL model to project cells and genes.
#' @returns
#' Data frame with plotting coordinates and basic information about type and
#' name.
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

#' Gets angle between plotting two directions.
#' Returns NA otherwise.
#' @param directions Row-wise matrix of directions.
#' @returns
#' Angle in degrees or NA if != 2 directions.
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

#' Flips APL direction if pointing away from majority of points.
#' @param apl_dirs Row-wise matrix of APL directions.
#' @param model The APL model to project cells and genes.
#' @param coords_type `prin_coords_cols` or `std_coords_cols`.
#' @inheritParams cluster_apl
#' @returns
#' The (flipped) APL directions.
.toggle_dir <- function(
    apl_dirs,
    caobj,
    cadir,
    model,
    coords_type = "prin_coords_cols"
) {
    # If the line points into the opposite direction of points
    # we flip the line.
    dapl_nms <- rownames(apl_dirs)
    ca_coords <- methods::slot(caobj, name = coords_type)
    for (d in seq_len(nrow(apl_dirs))) {
        sel <- match(
            names(cadir@cell_clusters)[cadir@cell_clusters == dapl_nms[d]],
            rownames(ca_coords)
        )

        if (length(sel) == 0) next
        cell_coords <- model(ca_coords)[sel, ]
        # cell_coords <- model(caobj@prin_coords_cols)[sel, ]

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

#' Adds information on the clusters to the plotting data frame.
#' @inheritParams cluster_apl
#' @param plot_df The plotting data frame.
#' @returns
#' The enhanced plotting data frame.
.add_cluster_info <- function(
    plot_df,
    cadir,
    cluster = NULL,
    highlight_cluster = FALSE,
    label_genes = FALSE) {
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
    sel_cells <- stats::na.omit(sel_cells)
    sel_genes <- match(names(cadir@gene_clusters)[gene_cl], plot_df$sample)
    sel_genes <- stats::na.omit(sel_genes)


    if (isTRUE(highlight_cluster)) {
        sel <- c() # empty in case we only want to plot lines.
        sel <- c(sel, sel_cells, sel_genes)
        sel <- stats::na.omit(sel)

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

        if (show_genes && label_genes) {
            plot_df$gene_score <- -Inf
            rnks <- cadir@gene_ranks[[cluster]]
            idx <- which(
                plot_df$sample %in% rnks$Rowname
            )
            ord <- base::match(plot_df[idx,]$sample, rnks$Rowname)
            rnks <- rnks[ord,]
            stopifnot(identical(rnks$Rowname, plot_df$sample[idx]))

            plot_df$gene_score[idx] <- rnks$Score
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

#' Creates ggplot for a cluster association plot
#' @param ggplt A ggplot2 object. It should be the result of ggplot() with
#' aes() defined but without any geoms_*.
#' @param plot_df The plotting data frame.
#' @inheritParams cluster_apl
#' @returns
#' A ggplot object. The basic cluster APL plot.
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

        ggplt <- .add_scales(
            ggplt = ggplt,
            point_size = point_size,
            size_factor = size_factor
        )

        if (isTRUE(highlight_cluster)) {
            ggplt <- .highlight_cluster(
                ggplt = ggplt,
                both = (show_cells && show_genes)
            )

            if (isTRUE(label_genes) && isTRUE(show_genes)) {
                ggplt <- .label_genes(
                    ggplt = ggplt,
                    plot_df = plot_df,
                    ntop = ntop,
                    show_cells = show_cells
                )
            }
        } else {
            ncls <- length(levels(plot_df$cluster))
            if (ncls <= 7) {
                ggplt <- ggplt +
                    scale_color_mpimg(name = "mpimg")
            } else if ( ncls > 7 && ncls <= 13) {
                ggplt <- ggplt +
                    scale_color_mpimg(name = "mpi_extend")
            }
        }
    }

    return(ggplt)
}

#' Adds size and shape scales to cluster APL ggplot.
#' @inheritParams .cluster_apl_points
#' @returns
#' ggplot with with size and shape scale added.
.add_scales <- function(ggplt, point_size = 1.5, size_factor = 2) {
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

#' Adds specific color schemes to highlight clusters.
#' @param ggplt A ggplot with geoms etc. already added.
#' @param both If TRUE highlights the cluster for both cells and genes,
#' otherwise cells only.
#' @returns
#' ggplot with scale_color added.
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

#' Helper function to label the top `ntop` genes of a cluster in the APL.
#' @inheritParams .cluster_apl_points
#' @param show_cells TRUE/FALSE whether to show cells or not.
#' @returns
#' ggplot object with labels added.
.label_genes <- function(ggplt, plot_df, show_cells, ntop = 15) {
    if (isTRUE(show_cells)){
        to_highlight <- (plot_df$cluster == "gene_cluster")
    } else {
        to_highlight <- (plot_df$cluster == "cluster")
    }

    dfh <- plot_df[to_highlight, ]
    dfh <- utils::head(dfh[order(dfh$gene_score, decreasing = TRUE), ], ntop)
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

#' Adds lines of other clusters to the APL.
#' @inheritParams cluster_apl
#' @param ggplt More or less finished ggplot object.
#' @param apl_dir The cluster directions in APL coordinates.
#' @param show_both TRUE if both cells and genes are shown.
#' @returns
#' ggplot with lines added.
.add_lines <- function(ggplt, show_lines, apl_dir, highlight_cluster = FALSE, show_both = FALSE) {
    if (is.character(show_lines)) {
        clusters <- show_lines
    } else {
        clusters <- rownames(apl_dir)
    }

    ltypes <- vector(mode = "character", length = length(clusters))
    names(ltypes) <- clusters
    slopes <- vector(mode = "numeric", length = length(clusters))
    names(slopes) <- clusters
    cl_name <- vector(mode = "character", length = length(clusters))
    names(cl_name) <- clusters

    for (c in clusters) {
        idx <- which(rownames(apl_dir) == c)
        is_x <- is_xaxis(apl_dir[idx, ])

        if (is_x) {
            ltype <- "solid"
            if (isTRUE(highlight_cluster)) {
                if (isTRUE(show_both)) {
                    hcluster <- "cell_cluster"
                } else {
                    hcluster <- "cluster"
                }
            } else {
                hcluster <- c
            }

        } else {
            ltype <- "dashed"
            if (isTRUE(highlight_cluster)) {
                if (isTRUE(show_both)) {
                    hcluster <- "cell_other"
                } else {
                    hcluster <- "other"
                }
            } else {
                hcluster  <- c
            }
        }

        i <- which(clusters == c)
        ltypes[i] <- ltype
        slopes[i] <- slope(lines = apl_dir[idx, ], dims = 1:2)
        cl_name[i] <- hcluster
    }

    ldf <- data.frame(
        "ltype" = as.factor(ltypes),
        "slopes" = slopes,
        "cluster" = as.factor(cl_name)
    )

    ldf$intercept <- 0
    ldf$orig_x <- 0
    ldf$orig_y <- 0

    ggplt <- ggplt + ggplot2::geom_abline(
        data = ldf,
        ggplot2::aes(
            slope = slopes,
            linetype = ltype,
            color = cluster,
            intercept = intercept
        ),
        linewidth = 1
    ) +
        ggplot2::scale_linetype_manual(
            values = c("dashed" = "dashed", "solid" = "solid")
        ) +
        ggplot2::geom_point(
            data = ldf,
            ggplot2::aes(
                x = orig_x,
                y = orig_y,
                color = cluster,
                text = NULL
            )
        )

    return(ggplt)
}

#' Checks if a direction is the x-axis (is the APL direction)
#' @param apl_dir A single direction.
#' @returns
#' TRUE if apl_dir is equivalent to the x-axis within tolerances,
#' FALSE otherwise.
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

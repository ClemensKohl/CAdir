#' Create a model to project points into an Association Plot
#' @param caobj A cacomp object.
#' @param direction Normed direction vector of the APL plot.
#' @param group A vector of indices which indicate the points
#' that belong to the cluster.
#' @returns
#' A model that can be used to project new points onto the APL plot.
apl_model <- function(
    caobj,
    direction,
    group = NULL
) {
    stopifnot(is(caobj, "cacomp"))

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
    }

    if (cosangle < 0) {
        avg_group_coords <- -avg_group_coords
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
#' @returns
#' An APL plot (ggplot2 object).
cluster_apl <- function(
        caobj,
        cadir,
        direction,
        group,
        cluster_id = "NA") {
    stopifnot(is(caobj, "cacomp"))
    stopifnot(is(cadir, "cadir"))

    # ensure that clusters and directions are coherent
    cadir <- rename_clusters(cadir)
    if (nrow(cadir@directions) > 1) {

        ang <- min(
            rad2deg(get_angle(cadir@directions[1, ], cadir@directions[2, ])),
            rad2deg(get_angle(-cadir@directions[1, ], cadir@directions[2, ]))
        )
    } else {
        ang = 0
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
    df$cluster <- 0
    sel <- match(names(cadir@cell_clusters), df$sample)
    df$cluster[sel] <- cadir@cell_clusters
    df$cluster <- as.factor(df$cluster)


    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = cluster)) +
        ggplot2::geom_point() +
        ggplot2::geom_point(
            data = data.frame(x = 0, y = 0),
            ggplot2::aes(x, y), color = "red"
        ) +
        ggplot2::ggtitle(paste0(
            "Cluster: ",
            cluster_id,
            ", CA-angle: ",
            round(ang, 2)
        )) +
        ggplot2::theme_bw()

    for (d in seq_len(nrow(dapl))) {
        p <- p + ggplot2::geom_abline(
            intercept = 0,
            slope = slope(lines = dapl[d, ], dims = 1:2),
            color = "red",
            linetype = "dashed",
            size = 1
        )
    }

    return(p)
}

# FIXME: WIP
plot_clusters <- function() {
    stop("Not implemented")
}

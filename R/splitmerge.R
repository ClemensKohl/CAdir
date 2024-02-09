# TODO: Combine with code for sa_splitmerge

# TODO: Add documentation
# FIXME: Update to work with caclust class
sub_cluster <- function(clustering, points, idx) {
    sel <- which(clustering@cell_clusters == idx)
    cl <- points[sel, ]
    res <- dirclust(
        points = cl,
        k = 2,
        epochs = 15,
        init = "kmeanspp",
        lines = NULL,
        log = FALSE
    )

    return(res)
}

# TODO: Add documentation
decide_split <- function(directions, cutoff = 30) {
    cutoff <- 1 - (deg2rad(cutoff) / pi)

    angsim <- max(
        get_ang_sim(directions[1, ], directions[2, ]),
        get_ang_sim(-directions[1, ], directions[2, ])
    )

    to_split <- angsim < cutoff

    return(to_split)
}
# FIXME: REMOVE OR FINISH.
add_split_cluster <- function(cadir, new_cl){

    # Update the clustering results and add the new directions.
    new_cl <- sort(unique(sres@cell_clusters))
    extra <- max(clustering@cell_clusters) + 1

    sel <- which(clustering@cell_clusters == i)
    cl1 <- sres@cell_clusters == new_cl[1]
    cl2 <- sres@cell_clusters == new_cl[2]

    if (!length(cl1) > 1 && length(cl2) > 1) next

    sres@cell_clusters[cl1] <- i
    sres@cell_clusters[cl2] <- extra

    clustering@cell_clusters[sel] <- sres@cell_clusters

    rownames(sres@directions) <- paste0("line", c(i, extra))
    clustering@directions[paste0("line", i), ] <- sres@directions[paste0("line", i), ]

    rn <- rownames(clustering@directions)
    clustering@directions <- rbind(
        clustering@directions,
        sres@directions[paste0("line", extra), ]
    )

    rownames(clustering@directions) <- c(rn, paste0("line", extra))

    return(cadir)
}

# TODO: Update names of clusters
# FIXME: Update to work with caclust class
# TODO: Find better way to deal with plots.
# TODO: Make function simpler.
determine_sub_clusters <- function(clustering,
    caobj,
    cutoff = 30,
    min_cells = 5,
    apl = FALSE,
    make_plots = FALSE) {
    cls <- sort(unique(clustering@cell_clusters))

    if (isTRUE(make_plots)) aplplots <- list()

    for (i in cls) {
        if (sum(clustering@cell_clusters == i) < 2) next

        sres <- sub_cluster(
            clustering = clustering,
            points = caobj@prin_coords_cols,
            idx = i
        )

        if (length(unique(sres@cell_clusters)) < 2) next
        stopifnot(length(unique(sres@cell_clusters)) == 2)


        elems <- as.numeric(table(sres@cell_clusters))

        if (any(elems < min_cells)) next

        dir <- sres@directions

        to_split <- decide_split(dir, cutoff = cutoff)

        if (isTRUE(to_split)) {
            message(paste0("Splitting cluster ", i))

            if (isTRUE(make_plots)) {
                p <- cluster_apl(
                    caobj = caobj,
                    cadir = sres,
                    apl_dir = clustering@directions[i, ],
                    indx_group = which(clustering@cell_clusters == i),
                    cluster_id = i
                )

                aplplots[[paste0("cluster_", i)]] <- p
            }

            # Update the clustering results and add the new directions.
            new_cl <- sort(unique(sres@cell_clusters))
            extra <- max(clustering@cell_clusters) + 1

            sel <- which(clustering@cell_clusters == i)
            cl1 <- sres@cell_clusters == new_cl[1]
            cl2 <- sres@cell_clusters == new_cl[2]

            if (!length(cl1) > 1 && length(cl2) > 1) next

            sres@cell_clusters[cl1] <- i
            sres@cell_clusters[cl2] <- extra

            clustering@cell_clusters[sel] <- sres@cell_clusters

            rownames(sres@directions) <- paste0("line", c(i, extra))
            clustering@directions[paste0("line", i), ] <- sres@directions[paste0("line", i), ]

            rn <- rownames(clustering@directions)
            clustering@directions <- rbind(
                clustering@directions,
                sres@directions[paste0("line", extra), ]
            )

            rownames(clustering@directions) <- c(rn, paste0("line", extra))
        }
    }

    if (isTRUE(make_plots)) {
        clustering$split_plots <- aplplots # FIXME: not a valid slot
    }

    clustering <- rename_clusters(clustering)

    return(clustering)
}

# FIXME: Update to work with caclust class
cameans_splitmerge <- function(caobj,
                               k,
                               cutoff = 40,
                               min_cells = 5,
                               epochs = 15,
                               reps = 5,
                               make_plots = FALSE) {
    out <- dirclust(
        points = caobj@prin_coords_cols,
        k = k,
        epochs = epochs,
        init = "kmeanspp",
        lines = NULL,
        log = FALSE
    )

    if (isTRUE(make_plots)) {
        split_plots <- list()
        merge_plots <- list()
    }


    for (i in seq_len(reps)) {
        message("Iteration ", i)

        subout <- determine_sub_clusters(
            clustering = out,
            caobj = caobj,
            cutoff = cutoff,
            min_cells = min_cells,
            apl = FALSE,
            make_plots = make_plots
        )

        mergeout <- merge_clusters(
            caobj = caobj,
            cameans = subout,
            cutoff = cutoff,
            make_plots = make_plots
        )

        # Final k-means to refine new clusters.

        out <- CAkmeans(caobj@prin_coords_cols,
            k = ncol(mergeout$direction),
            lines = mergeout$direction,
            epochs = 5,
            log = FALSE
        )

        if (isTRUE(make_plots)) {
            split_plots[[paste0("rep_", i)]] <- subout$split_plots

            merge_plots[[paste0("rep_", i)]] <- mergeout$merge_plots
        }
    }

    if (isTRUE(make_plots)) {
        out$split_plots <- split_plots
        out$merge_plots <- merge_plots
    }
    return(out)
}

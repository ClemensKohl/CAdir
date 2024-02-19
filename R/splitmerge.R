# TODO: Combine with code for sa_splitmerge

#' Split a cluster into 2 sub-clusters.
#' @inheritParams dirclust
#' @param idx The cluster index to split.
#' @param cadir A `cadir` object.
#' @returns
#' A `cadir` object that contains only the two new clusters.
sub_cluster <- function(cadir, points, idx) {

    sel <- which(cadir@cell_clusters == idx)
    cl <- points[sel, ]
    if (length(cl) == 0) {
        stop("No points in cluster")
    }
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


#' Determine if a cluster should be split.
#' @param directions A row vector matrix of directions.
#' @param cutoff The cutoff angle in degrees.
#' @returns
#' A logical indicating if the cluster should be split.
decide_split <- function(directions, cutoff = deg2rad(30)) {
    cutoff <- 1 - (cutoff / pi)

    angsim <- max(
        get_ang_sim(directions[1, ], directions[2, ]),
        get_ang_sim(-directions[1, ], directions[2, ])
    )

    to_split <- angsim < cutoff

    return(to_split)
}

# TODO: Make function smarter.
# Only check clusters you haven't checked yet.

# TODO: Fix documentation

#' Split clusters based on angle.
#' @inheritParams dirclust_splitmerge
#' @param cadir A `cadir` object.
#' @returns A `cadir` object with split clusters.
split_clusters <- function(
    cadir,
    caobj,
    cutoff = 30,
    method = "random",
    min_cells = 5,
    make_plots = FALSE,
    counts = NULL,
    reps = 100,
    apl_quant = 0.99
) {

    cls <- sort(unique(cadir@cell_clusters))

    for (i in cls) {

        if (sum(cadir@cell_clusters == i) < 2) next

        sres <- sub_cluster(
            cadir = cadir,
            points = caobj@prin_coords_cols,
            idx = i
        )

        if (length(unique(sres@cell_clusters)) < 2) next
        stopifnot(length(unique(sres@cell_clusters)) == 2)

        elems <- as.numeric(table(sres@cell_clusters))
        if (any(elems < min_cells)) next

        # If the cutoff is set to NULL we computer the S-alpha cutoff.
        if (is.null(cutoff)) {

            grp_idx <- which(cadir@cell_clusters == i)

            aplcds <- apl_dir_coords(
                cadir = sres,
                caobj = caobj,
                apl_dir = cadir@directions[f2n(i), ],
                group = grp_idx
            )

            if (method == "permutation" && is.null(counts)) {
                warning(
                    "No count matrix for permutation supplied.",
                    "Switching to random directions method."
                )
                counts <- NULL
                method <- "random"
            }

            cutoff <- get_apl_cutoff(
                caobj = caobj,
                counts = counts,
                method = method,
                apl_cols = aplcds$apl_cols,
                group = grp_idx,
                dims = caobj@dims,
                quant = apl_quant
            )

            to_split <- decide_split(aplcds$apl_dirs, cutoff = cutoff)

        } else {
            to_split <- decide_split(sres@directions, cutoff = cutoff)
        }


        if (isTRUE(to_split)) {
            message(paste0("Splitting cluster ", i))

            if (isTRUE(make_plots)) {

                p <- cluster_apl(
                    caobj = caobj,
                    cadir = sres,
                    direction = cadir@directions[i, ],
                    group = which(cadir@cell_clusters == i),
                    cluster_id = i
                )
                rep <- paste0("rep_", cadir@log$last_rep)
                nms_rep <- names(cadir@plots$splits[[rep]])
                nm <- paste("cluster", i, collapse = "_")

                if (nm %in% names(nms_rep)) {
                    cnt <- sum(grepl(nm, nms_rep))
                    nm <- paste0(nm, cnt + 1)
                }

                cadir@plots$splits[[rep]][[nm]] <- p
            }

            # Update the clustering results and add the new directions.
            new_cl <- sort(unique(sres@cell_clusters))
            lvls <- as.numeric(as.character(levels(cadir@cell_clusters)))
            extra <- max(lvls) + 1

            sel <- which(cadir@cell_clusters == factor(i, levels = lvls))
            cl1 <- sres@cell_clusters == new_cl[1]
            cl2 <- sres@cell_clusters == new_cl[2]

            if (!length(cl1) > 1 && length(cl2) > 1) next

            new_lvls <- c(lvls, extra)

            # update the factor levels
            levels(sres@cell_clusters) <- new_lvls
            sres@cell_clusters[cl1] <- i
            sres@cell_clusters[cl2] <- extra

            levels(cadir@cell_clusters) <- new_lvls
            cadir@cell_clusters[sel] <- sres@cell_clusters

            rownames(sres@directions) <- paste0("line", c(i, extra))
            cadir@directions[paste0("line", i), ] <-
                sres@directions[paste0("line", i), ]

            rn <- rownames(cadir@directions)
            cadir@directions <- rbind(
                cadir@directions,
                sres@directions[paste0("line", extra), ]
            )

            rownames(cadir@directions) <- c(rn, paste0("line", extra))

            cadir <- rename_clusters(cadir)

            cadir <- split_clusters(
                cadir = cadir,
                caobj = caobj,
                cutoff = cutoff,
                min_cells = min_cells,
                make_plots = make_plots
            )
            return(cadir)

        }
    }
    return(cadir)
}

# TODO: handle plots better and make recursive smoother.
# TODO: replace tmp variables with class slots if possible.
# TODO: Fix documentation

#' Merge clusters based on angle
#' @inheritParams dirclust_splitmerge
#' @param cadir A `cadir` object.
#' @returns A `cadir` object with merged clusters.
merge_clusters <- function(caobj,
                           cadir,
                           cutoff,
                           method = "random",
                           counts = NULL,
                           apl_quant = 0.99,
                           make_plots = FALSE) {

    samples <- caobj@prin_coords_cols
    clusters <- cadir@cell_clusters
    directions <- cadir@directions

    if (isTRUE(make_plots)) {
        rep <- paste0("rep_", cadir@log$last_rep)

        if (length(cadir@plots$merges[[rep]]) == 0) {
            mergeplots <- list()
        } else {
            mergeplots <- cadir@plots$merges[[rep]]
        }
    }

    if (is.null(cutoff)) {
        candidates <- get_apl_mergers(cadir = cadir,
                                      caobj = caobj,
                                      method = method,
                                      counts = counts,
                                      apl_quant = apl_quant)


        candidates <- apply(candidates, 1, function(x) which(x))

    } else {
        asim_cutoff <- 1 - cutoff / pi
        sim <- get_ang_sim(directions, directions)
        # Set the lower diagonal to 0
        sim[lower.tri(sim, diag = TRUE)] <- 0
        candidates <- apply(sim, 1, function(x) which(x >= asim_cutoff))
    }

    sel <- which(lengths(candidates) > 0)

    for (s in sel) {
        cds <- candidates[[s]]
        mergers <- c(s, cds)

        message(paste0("Merging cluster ", s, " with ", cds, "\n"))

        cls <- which(clusters %in% mergers)
        if (length(cls) == 0) next

        new_dir <- drop(total_least_squares(samples[cls, , drop = FALSE]))

        if (isTRUE(make_plots)) {
            sres <- new("cadir",
                cell_clusters = clusters[cls],
                directions = directions[mergers, ],
            )

            sres@cell_clusters <- match(
                sres@cell_clusters,
                sort(unique(sres@cell_clusters))
            )
            names(sres@cell_clusters) <- names(clusters[cls])

            p <- cluster_apl(
                caobj = caobj,
                cadir = sres,
                direction = new_dir,
                group = cls,
                cluster_id = s
            )

            nm <- paste("cluster", mergers, collapse = "_")

            if (nm %in% names(mergeplots)) nm <- paste0(nm, "_v2")
            mergeplots[[nm]] <- p
        }

        clusters[cls] <- s
        directions[s, ] <- new_dir
        directions <- directions[-cds, ]

        cadir@cell_clusters <- clusters
        cadir@directions <- directions

        if (length(sel) > 1) {
            cadir <- rename_clusters(cadir = cadir)

            if (isTRUE(make_plots)) {
                rep <- paste0("rep_", cadir@log$last_rep)
                cadir@plots$mergers[[rep]] <- append(
                    cadir@plots$mergers[[rep]],
                    mergeplots
                )
            }

            out <- merge_clusters(
                caobj = caobj,
                cadir = cadir,
                cutoff = cutoff,
                method = method,
                counts = counts,
                apl_quant = apl_quant,
                make_plots = make_plots
            )
            return(out)
        }
    }

    cadir@directions <- directions
    cadir@cell_clusters <- clusters
    cadir <- rename_clusters(cadir)

    if (isTRUE(make_plots)) {
        rep <- paste0("rep_", cadir@log$last_rep)
        cadir@plots$mergers[[rep]] <- append(
            cadir@plots$mergers[[rep]],
            mergeplots
        )
    }

    return(cadir)
}

# TODO: Fix documentation

#' Perform Clustering by CA directions with splitting and merging.
#' @inheritParams dirclust
#' @param caobj A `caclust` object.
#' @param reps Number of repetitions to perform the splitting and merging.
#' @param min_cells Minimum number of cells to form a cluster.
#' @param make_plots Logical. If `TRUE` plots are generated for each
#' split and merge
#' @param cutoff Degrees. The cutoff angle to split and merge clusters.
#' @param qcutoff The quantile cutoff for gene selection.
#' @return A `cadir` object with cell clusters.
dirclust_splitmerge <- function(caobj,
                                k,
                                cutoff = 40,
                                qcutoff = 0.8,
                                method = "random",
                                counts = NULL,
                                apl_quant = 0.99,
                                min_cells = 5,
                                epochs = 15,
                                reps = 5,
                                make_plots = FALSE) {
    # Convert cutoff to radians
    if (!is.null(cutoff)) cutoff <- deg2rad(cutoff)

    out <- dirclust(
        points = caobj@prin_coords_cols,
        k = k,
        epochs = epochs,
        init = "kmeanspp",
        lines = NULL,
        log = FALSE
    )

    for (i in seq_len(reps)) {
        message("Iteration ", i)

        out@log$last_rep <- i

        out <- split_clusters(
            cadir = out,
            caobj = caobj,
            cutoff = cutoff,
            method = method,
            counts = counts,
            apl_quant = apl_quant,
            min_cells = min_cells,
            make_plots = make_plots
        )

        out <- merge_clusters(
            caobj = caobj,
            cadir = out,
            cutoff = cutoff,
            method = method,
            counts = counts,
            apl_quant = apl_quant,
            make_plots = make_plots
        )

        # Final k-means to refine new clusters.

        out <- dirclust(
            points = caobj@prin_coords_cols,
            k = ncol(out@directions),
            lines = out@directions,
            epochs = 5,
            log = FALSE
        )
    }

    out@gene_clusters <- assign_genes(
        caobj = caobj,
        cadir = out,
        qcutoff = qcutoff,
        coords = "prin"
    )

    out <- rename_clusters(out)
    return(out)
}

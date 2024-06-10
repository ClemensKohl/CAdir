# TODO: make a verbose toggle.

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
    apl_cutoff_reps = 100,
    apl_quant = 0.99
) {

    fun_args <- match.call()
    cls <- sort(unique(cadir@cell_clusters))

    for (i in cls) {
        # cat("Type of i", class(i))
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

            cutoff_exists <- is_stored(cadir = cadir, fun_args = fun_args)

            if (isTRUE(cutoff_exists)) {
                cutoff <- cadir@parameters$sa_cutoff
            } else {
                cutoff <- get_apl_cutoff(
                    caobj = caobj,
                    counts = counts,
                    method = method,
                    group = grp_idx,
                    quant = apl_quant,
                    apl_cutoff_reps = apl_cutoff_reps
                )
                cadir@parameters$sa_cutoff <- cutoff
            }

            to_split <- decide_split(aplcds$apl_dirs, cutoff = cutoff)

        } else {
            to_split <- decide_split(sres@directions, cutoff = cutoff)
        }


        if (isTRUE(to_split)) {
            message(paste0("\tSplitting cluster ", i))

            if (isTRUE(make_plots)) {

                p <- cluster_apl(
                    caobj = caobj,
                    cadir = sres,
                    direction = cadir@directions[f2n(i), ],
                    group = which(cadir@cell_clusters == i),
                    cluster = NULL
                )

                rep <- paste0("rep_", cadir@log$last_rep)
                nms_rep <- names(cadir@plots$splits[[rep]])
                nm <- paste("cluster", f2n(i), collapse = "_", sep = "")

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

            cadir@distances <- matrix(0, 0, 0)

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
                           make_plots = FALSE,
                           apl_cutoff_reps = 100) {
    samples <- caobj@prin_coords_cols
    clusters <- cadir@cell_clusters
    directions <- cadir@directions

    if (is.null(cutoff)) {
        candidates <- get_apl_mergers(cadir = cadir,
                                      caobj = caobj,
                                      method = method,
                                      counts = counts,
                                      apl_quant = apl_quant,
                                      apl_cutoff_reps = apl_cutoff_reps)


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

        message(paste0("\tMerging cluster ", s, " with ", cds))

        cls <- which(clusters %in% mergers)
        if (length(cls) == 0) next

        new_dir <- drop(total_least_squares(samples[cls, , drop = FALSE]))

        if (isTRUE(make_plots)) {

            rep <- paste0("rep_", cadir@log$last_rep)

            sres <- methods::new("cadir",
                cell_clusters = clusters[cls],
                directions = directions[mergers, ]
            )

            names(sres@cell_clusters) <- names(clusters[cls])

            p <- cluster_apl(
                caobj = caobj,
                cadir = sres,
                direction = new_dir,
                group = cls,
                cluster = s
            )

            nms_rep <- names(cadir@plots$merges[[rep]])
            nm <- paste("cluster", mergers, collapse = "_", sep = "")

            if (nm %in% names(nms_rep)) {
                cnt <- sum(grepl(nm, nms_rep))
                nm <- paste0(nm, cnt + 1)
            }

            cadir@plots$merges[[rep]][[nm]] <- p
        }

        clusters[cls] <- s
        directions[s, ] <- new_dir
        directions <- directions[-cds, , drop = FALSE]


        cadir@cell_clusters <- clusters
        cadir@directions <- directions
        cadir@distances <- matrix(0, 0, 0) # dists not true anymore.
        cadir@cl2dir[[cl2nm(s)]] <- s # ensure dict is up to date.

        cadir <- rename_clusters(cadir = cadir)

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

    cadir@directions <- directions
    cadir@cell_clusters <- clusters
    cadir <- rename_clusters(cadir)

    return(cadir)
}


#' Perform Clustering by CA directions with splitting and merging.
#' @inheritParams dirclust
#' @inheritParams get_apl_cutoff
#' @param caobj A `caclust` object.
#' @param reps Number of repetitions to perform the splitting and merging.
#' @param min_cells Minimum number of cells to form a cluster.
#' @param make_plots Logical. If `TRUE` plots are generated for each
#' @param apl_quant The quantile to use for the APL cutoff. Only used
#' if cutoff is `NULL`.
#' split and merge
#' @param cutoff Degrees. The cutoff angle to split and merge clusters.
#'  If `NULL` the cutoff angle is calculated based on the cutoff angle based
#'  on the quantile defined through `apl_quant`.
#' @param qcutoff The quantile cutoff for gene selection.
#' @return A `cadir` object with cell clusters.
#' @seealso [get_apl_cutoff()]
#' @export
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
                                apl_cutoff_reps = 100,
                                make_plots = FALSE) {
    # Convert cutoff to radians
    if (!is.null(cutoff)) cutoff <- deg2rad(cutoff)
    fun_args <- match.call()

    # Set up logging.
    log <- list()
    log$clusters <- as.data.frame(matrix(0,
                                   ncol = 1,
                                   nrow = nrow(caobj@prin_coords_cols)))
    colnames(log$clusters) <- "root"

    log$directions <- as.data.frame(matrix(0,
                                    ncol = ncol(caobj@prin_coords_cols),
                                    nrow = 1))
    colnames(log$directions) <- colnames(caobj@prin_coords_cols)
    log$directions <- cbind(iter = "root",
                            dirname = "root",
                            log$directions)

    # Initial CAdir clustering.
    out <- dirclust(
        points = caobj@prin_coords_cols,
        k = k,
        epochs = epochs,
        init = "kmeanspp",
        lines = NULL,
        log = FALSE
    )

    out <- rename_clusters(out)

    log <- log_iter(log = log,
                    cadir = out,
                    name = "iter_0")

    for (i in seq_len(reps)) {
        message("Iteration ", i)

        iter_nm <- paste0("iter_", i)
        out@log$last_rep <- i

        out <- split_clusters(
            cadir = out,
            caobj = caobj,
            cutoff = cutoff,
            method = method,
            counts = counts,
            apl_quant = apl_quant,
            min_cells = min_cells,
            make_plots = make_plots,
            apl_cutoff_reps = apl_cutoff_reps
        )

        log <- log_iter(log = log,
                        cadir = out,
                        name = paste0("split_", iter_nm))

        out <- dirclust(
            points = caobj@prin_coords_cols,
            k = ncol(out@directions),
            epochs = 5,
            log = FALSE,
            cadir = out
        )

        out <- rename_clusters(out)


        log <- log_iter(log = log,
                        cadir = out,
                        name = paste0("interS_", iter_nm))

        out <- merge_clusters(
            caobj = caobj,
            cadir = out,
            cutoff = cutoff,
            method = method,
            counts = counts,
            apl_quant = apl_quant,
            make_plots = make_plots,
            apl_cutoff_reps = apl_cutoff_reps
        )

        log <- log_iter(log = log,
                        cadir = out,
                        name = paste0("merge_", iter_nm))

        out <- dirclust(
            points = caobj@prin_coords_cols,
            k = ncol(out@directions),
            epochs = 5,
            log = FALSE,
            cadir = out
        )
        out <- rename_clusters(out)

        log <- log_iter(log = log,
                        cadir = out,
                        name = paste0("interM_", iter_nm))

    }

    out@gene_clusters <- assign_genes(
        caobj = caobj,
        cadir = out,
        qcutoff = qcutoff,
        coords = "prin"
    )

    out <- rename_clusters(out)

    log <- log_iter(log = log,
                    cadir = out,
                    name = paste0("end", iter_nm))


    to_keep <- (!names(out$log) %in% names(log))
    out@log <- c(out@log[to_keep], log)
    out@parameters$qcutoff <- qcutoff
    out@parameters$call <- fun_args


    return(out)
}

#' @include classes.R
#' @importFrom CAbiNet annotate_biclustering
#' @importFrom CAbiNet annotate_by_goa
NULL

# FIXME: Add documentation
setMethod(
    f = "annotate_by_goa",
    signature = (object <- "cadir"),
    function(obj,
             goa_res,
             alpha = 0.05) {
        stopifnot(is(obj, "cadir"))

        # cell clusters
        ccs <- cell_clusters(obj)
        ccs_nm <- names(ccs)

        unccs <- sort(unique(ccs))

        ccs <- as.character(ccs)
        unccs <- as.character(unccs)

        # gene clusters
        gcs <- gene_clusters(obj)
        gcs_nm <- names(gcs)
        ungcs <- sort(unique(gcs))

        gcs <- as.character(gcs)
        ungcs <- as.character(ungcs)

        # Combine clusters
        allcs <- sort(unique(c(unccs, ungcs)))

        # Directions
        dirs <- obj@directions
        stopifnot(nrow(dirs) == length(allcs))

        # Create dataframe with all goa results
        nclust <- length(goa_res)

        goa_res <- lapply(goa_res, function(x, nc = nclust) {
            subs <- min(nrow(x), nc)
            x[seq_len(subs), ]
        })

        # Solve assignment problem with the hungarian algorithm.
        cluster_anno <- CAbiNet::assign_cts(goa_res)

        # Rename clusters based on GSE.
        for (c in seq_len(length(allcs))) {
            noanno <- paste0("cluster", allcs[c])

            if (!allcs[c] %in% cluster_anno$cluster) {
                ct <- noanno
            } else {
                ct <- cluster_anno[cluster_anno$cluster == allcs[c], "cell_type"]
                padj <- cluster_anno[cluster_anno$cluster == allcs[c], "padj"]
                if (padj > alpha) ct <- noanno
            }

            if (allcs[c] %in% ungcs) {
                sel <- which(gcs == allcs[c])

                if (length(ct) == 0) {
                    gcs[sel] <- noanno
                } else {
                    gcs[sel] <- ct
                }
            }

            if (allcs[c] %in% unccs) {
                sel <- which(ccs == allcs[c])

                if (length(ct) == 0) {
                    ccs[sel] <- noanno
                } else {
                    ccs[sel] <- ct
                }
            }

            if (c <= nrow(dirs)) {
                if (length(ct) == 0) {
                    rownames(dirs)[c] <- noanno
                } else {
                    rownames(dirs)[c] <- ct
                }
            }
        }

        names(gcs) <- gcs_nm
        names(ccs) <- ccs_nm

        # update results
        obj@cell_clusters <- as.factor(ccs)
        obj@gene_clusters <- as.factor(gcs)
        obj@directions <- dirs

        stopifnot(validObject(obj))
        return(obj)
    }
)


# TODO: Add documentation
#' Annotate the biclustering
#' @rdname annotate_biclustering
#' @export
setMethod(
    f = "annotate_biclustering",
    signature = (obj <- "cadir"),
    function(obj,
             universe,
             org,
             set = "CellMarker",
             alpha = 0.05,
             min_size = 10,
             max_size = 500,
             ...) {
        stopifnot(is(obj, "caclust"))

        goa_res <- CAbiNet::per_cluster_goa(
            cabic = obj,
            universe = universe,
            set = set,
            org = org,
            min_size = min_size,
            max_size = max_size
        )

        cabic <- annotate_by_goa(
            obj = obj,
            goa_res = goa_res,
            alpha = alpha
        )

        return(cabic)
    }
)

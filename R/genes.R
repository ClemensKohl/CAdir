#' Returns the indices of all points with a norm outside of a sphere
#'  with radius of the defined quantile of the vector norm.
#'
#' @param x matrix of row vectors
#' @param qcutoff quantile.
#'
#' @returns
#' Indices of points lying outside of sphere
ca_sphere_idx <- function(x, qcutoff = 0.8) {
    xn <- row_norm(x)
    q <- stats::quantile(xn, qcutoff)
    idx <- which(xn > q)

    return(idx)
}

#' Assigns genes to cell clusters based on distance to closest line.
#' @description
#' `assign_genes` removes all genes inside of a sphere which's radius
#' is defined by the quantile cutoff on the norm of the genes (using
#'  the principal coordinates by default).
#' The remaining genes are then assigned to the closest line using the
#' standard coordinates.
#' @param caobj A cacomp object.
#' @param cadir A cadir object.
#' @param qcutoff The quantile cutoff for gene selection.
#' @param coords The coordinates to use for gene selection ("std" or "prin").
#' @returns
#' Gene clusters.
#' @export
assign_genes <- function(caobj,
                         cadir,
                         qcutoff = NULL,
                         coords = "prin") {
    if (coords == "prin") {
        idx <- ca_sphere_idx(caobj@prin_coords_rows, qcutoff = qcutoff)
    } else if (coords == "std") {
        idx <- ca_sphere_idx(caobj@std_coords_rows, qcutoff = qcutoff)
    } else {
        stop("Invalid coords argument")
    }

    X <- caobj@std_coords_rows[idx, ]
    gene_nms <- rownames(X)

    ldist <- dist_to_line(X, cadir@directions, row_norm(X))
    # find closest line
    clusters <- apply(ldist, 1, which.min)

    cell_lvls <- as.numeric(as.character(levels(cadir@cell_clusters)))
    lvls <- sort(unique(c(unique(clusters), cell_lvls)))

    clusters <- factor(clusters, levels = lvls)
    names(clusters) <- gene_nms

    return(clusters)
}

# TODO: Add documentation
#' @export
rank_genes <- function(cadir, caobj) {

    # Step1: Get cutoff
    alpha <- cadir@parameters$sa_cutoff

    # Step 2: Get APL coordinates of co-clustered genes
    # for the cluster direction.
    gcs <- sort(unique(cadir@gene_clusters))
    gcs_lvls <- as.numeric(gcs)

    all_ranks <- list()
    for (c in gcs_lvls) {

        if (is.null(alpha)) {
            alpha <- get_apl_cutoff(
                caobj = caobj,
                method = "random",
                group = which(cadir@cell_clusters == gcs[c]),
                quant = 0.99,
                apl_cutoff_reps = 100
            )
        }

        direction <- cadir@directions[c, ]

        model <- apl_model(
            caobj = caobj,
            direction = direction,
            group = which(cadir@cell_clusters == gcs[c])
        )

        gene_coords <- model(caobj@prin_coords_rows)

        # subset genes to cluster genes
        gene_coords <- gene_coords[
            which(cadir@gene_clusters == gcs[c]), ,
            drop = FALSE
        ]

        # Step 3: Score by cutoff
        score <- gene_coords[, 1] - (gene_coords[, 2] * alpha)

        ranking <- data.frame(
            "Rowname" = rownames(gene_coords),
            "Score" = score,
            "Row_num" = seq_len(nrow(gene_coords)),
            "cluster" = gcs[c]
        )

        ranking <- ranking[order(ranking$Score, decreasing = TRUE), ]
        ranking$Rank <- seq_len(nrow(ranking))

        all_ranks[[as.character(gcs[c])]] <- ranking
        # all_ranks[[paste0("cluster_", c)]] <- ranking
    }
    # Step 4: Store results is new slot (needs to be created)
    cadir@gene_ranks <- all_ranks

    return(cadir)
}

#TODO: Add documentation
#' @export
top_genes <- function(cadir, cutoff = 0) {
    top_list <- list()
    gene_ranks <- cadir@gene_ranks
    for (i in seq_len(length(gene_ranks))) {
        relevant_genes <- gene_ranks[[i]][gene_ranks[[i]]$Score > cutoff, ]
        del <- which(colnames(relevant_genes) == "Rank")
        top_list[[names(gene_ranks)[i]]] <- relevant_genes[, -del]
    }

    topdf <- do.call("rbind", top_list)
    return(topdf)
}

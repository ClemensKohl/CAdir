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

#TODO: Implement per-cluster gene ranking by S-alpha score

rank_genes <- function(cadir, caobj){
    stop("Not implemented yet")

    # Step1: Get cutoff
    alpha <- attr(caobj@permuted_data, "cutoff")

    # Step 2: Get APL coordinates of co-clustered genes for the cluster direction.
    # Step 3: Score by cutoff
    score <- caobj@apl_rows[,1] - (caobj@apl_rows[,2] * cutoff_cotan)
    ranking <- data.frame("Rowname" = rownames(caobj@apl_rows),
                          "Score" = score,
                          "Row_num" = seq_len(nrow(caobj@apl_rows)))

    ranking <- ranking[order(ranking$Score, decreasing = TRUE),]
    ranking$Rank <- seq_len(nrow(ranking))

    # Step 4: Store results is new slot (needs to be created)

    return(cadir)
}

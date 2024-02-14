# TODO: Add documentation
ca_sphere_idx <- function(x, qcutoff = 0.8) {

    xn <- row_norm(x)
    q <- quantile(xn, qcutoff)
    idx <- which(xn > q)

    return(idx)
}

# TODO: Add documentation
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

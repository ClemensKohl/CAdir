one_to_rule_them_all <- function(caobj,
                                 cadir,
                                 cluster = NULL,
                                 direction,
                                 group,
                                 show_cells = TRUE,
                                 show_genes = FALSE,
                                 show_lines = TRUE,
                                 highlight_cluster = TRUE,
                                 # colour_by_group = FALSE,
                                 label_genes = FALSE,
                                 point_size = 1.5,
                                 size_factor = 2,
                                 ntop = 15) {
    stopifnot(methods::is(caobj, "cacomp"))
    stopifnot(methods::is(cadir, "cadir"))


    all_cls <- unique(c(
        levels(cadir@cell_clusters),
        levels(cadir@gene_clusters)
    ))

    if (!is.null(cluster)) {
        cluster <- as.character(cluster)
        if (!cluster %in% all_cls) cluster <- NULL
    }

    if (is.null(cluster)) {
        # Kinda redundant. placeholder if I want to do deal with special case.
        cell_grp <- seq_len(length(cadir@cell_clusters))
        gene_grp <- seq_len(length(cadir@gene_clusters))
    } else {
        cell_grp <- which(cadir@cell_clusters == cluster)
        gene_grp <- which(cadir@gene_clusters == cluster)
    }

    # Calculate angle if only two directions, NA otherwise
    ang <- .get_plot_angle(directions = cadir@directions)

    model <- apl_model(
        caobj = caobj,
        direction = direction,
        group = group
    )

    dapl <- model(cadir@directions)
    dapl_nms <- rownames(dapl)

    bool_sum <- show_cells + show_genes + show_cells
    coords <- list(NULL, "prin_coords_cols", "std_coords_cols")

    df <- .construct_df(
        cells = slot(caobj, name = coords[[bool_sum]]),
        genes = ifelse(show_genes, yes = caobj@prin_coords_rows, no = NULL),
        model = model
    )
}



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

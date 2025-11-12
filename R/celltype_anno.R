#' @include classes.R
#' @importFrom CAbiNet annotate_biclustering
#' @importFrom CAbiNet annotate_by_goa
NULL

#' Assign cell types to clusters using the Hungarian algorithm.
#' @description
#' Uses the hungarian algorithm (assignment problem)
#' to assign a cell type from the gene set overrepresentation
#' analysis to one (and only one) cluster.
#' @param goa_res List of goa results for each bicluster.
#' @returns
#' A data frame with the assigned cell types and adjusted p-values.
#' @export
assign_cts_logpval <- function(goa_res) {
  # Solve assignment problem with the hungarian algorithm.
  # Build cost matrix.
  goa_res <- dplyr::bind_rows(goa_res, .id = "cluster")

  cost_mat <- stats::reshape(
    data = goa_res[, c("cluster", "gene_set", "padj")],
    direction = "wide",
    idvar = "cluster",
    timevar = "gene_set",
    new.row.names = seq_len(length(unique(goa_res$cluster)))
  )

  cost_mat[is.na(cost_mat)] <- 1
  colnames(cost_mat) <- gsub("padj.", "", colnames(cost_mat), fixed = TRUE)

  if ("No_Cell_Type_Found" %in% colnames(cost_mat)) {
    rm_col <- which(colnames(cost_mat) == "No_Cell_Type_Found")
    cost_mat <- cost_mat[, -rm_col, drop = FALSE]
  }

  if (ncol(cost_mat) == 1) {
    stop(
      "GOA results do not contain any cell types! Check if any genes are in the gene sets!"
    )
  }

  clusters <- as.character(cost_mat$cluster)
  cell_types <- colnames(cost_mat)[2:ncol(cost_mat)]

  cost_mat <- as.matrix(cost_mat[, 2:ncol(cost_mat)], drop = FALSE)

  # solve assignment problem
  assignments <- RcppHungarian::HungarianSolver(cost_mat)$pairs

  assignments <- assignments[assignments[, 2] > 0, ]

  cluster_anno <- data.frame(
    cluster = clusters[assignments[, 1]],
    cell_type = cell_types[assignments[, 2]],
    padj = cost_mat[assignments]
  )

  return(cluster_anno)
}

#' Annotate CAbiNet results by gene overrepresentation
#'  analysis results.
#'
#' @description
#' `annotate_by_goa` takes a biclustering results such as outputted by `caclust`
#' and annotates it with the gene overrepresentation analysis results (goa).
#'
#' @inheritParams CAbiNet::annotate_by_goa
#' @param obj `cadir` object with biclustering results. Alternatively could be
#' a caclust object.
#'
#' @description
#' Conflicts between clusters that have the
#'  same highest ranking gene set are solved
#'  with the Hungarian/Munkres algorithm.
#'
#' @returns
#' Object of the same type as `obj`.
#'
#' @export
setMethod(
  f = "annotate_by_goa",
  signature = (object <- "cadir"),
  function(obj, goa_res, alpha = 0.05) {
    stopifnot(methods::is(obj, "cadir"))

    # dictionary
    dict <- obj@dict

    # cell clusters
    ccs <- cell_clusters(obj)
    ccs_nm <- names(ccs)
    unccs <- levels(ccs)

    # gene clusters
    gcs <- gene_clusters(obj)
    gcs_nm <- names(gcs)
    ungcs <- levels(gcs)

    # Combine clusters
    allcs <- unique(c(unccs, ungcs))

    # Convert to character
    ccs <- as.character(ccs)
    unccs <- as.character(unccs)
    gcs <- as.character(gcs)
    ungcs <- as.character(ungcs)
    allcs <- as.character(allcs)

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
    cluster_anno <- assign_cts_logpval(goa_res)

    # Rename clusters based on GSE.
    for (c in seq_len(length(allcs))) {
      noanno <- allcs[c]

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

          # Update name and entry in the dictionary
          dict_idx <- which(names(dict) == search_dict(dict, c))
          names(dict)[dict_idx] <- ct
          dict[[ct]] <- c
        }
      }
    }

    names(gcs) <- gcs_nm
    names(ccs) <- ccs_nm
    lvls <- rownames(dirs)

    # update results
    obj@cell_clusters <- factor(ccs, levels = lvls)
    obj@gene_clusters <- factor(gcs, levels = lvls)
    obj@directions <- dirs
    obj@dict <- dict

    stopifnot(methods::validObject(obj))
    return(obj)
  }
)

#' Annotate the biclustering
#' @inheritParams CAbiNet::annotate_biclustering
#' @param obj A cadir object, or,
#' alternatively a `caclust` or `SingleCellExperiment` object with biclustering
#' information.
#' @rdname annotate_biclustering
#' @export
setMethod(
  f = "annotate_biclustering",
  signature = (obj <- "cadir"),
  function(
    obj,
    universe,
    org,
    set = "CellMarker",
    alpha = 0.05,
    min_size = 10,
    max_size = 500,
    ...
  ) {
    stopifnot(methods::is(obj, "caclust"))

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

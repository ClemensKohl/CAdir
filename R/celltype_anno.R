#' @include classes.R
#' @importFrom CAbiNet annotate_biclustering
#' @importFrom CAbiNet annotate_by_goa
NULL

#' Filter gene sets by size
#' @description
#' This function cuts removes gene sets below and above a certain size.
#' @param gene_sets A named list of gene sets (name is the gene set,
#' each element of the list contains a vector of genes)
#' @param min_size Min. number of genes in the gene set.
#' @param max_size Max number of genes in the gene set.
#' Set to Inf if you want to keep all genes.
#' @returns
#' Filtered list of gene sets.
filter_gene_sets <- function(
  gene_sets,
  min_size = 10,
  max_size = 500,
  filter_literature = FALSE
) {
  if (is.na(min_size) || is.null(min_size)) {
    min_size <- 1
  }
  if (is.na(max_size) || is.null(max_size)) {
    max_size <- Inf
  } #.Machine$integer.max

  ## index of gene_sets in used.
  ## logical
  gene_sets_length <- lengths(gene_sets)
  idx <- min_size <= gene_sets_length & gene_sets_length <= max_size

  gene_sets <- gene_sets[idx]
  if (isTRUE(filter_literature)) {
    gene_sets <- gene_sets[!grepl("et_al\\.", names(gene_sets))]
  }
  return(gene_sets)
}

# NOTE: use fsgea package
perform_gsea <- function(
  salpha,
  gene_sets,
  min_size = 15,
  max_size = 500,
  set = "CellMarker"
) {
  gsea <- fgsea::fgsea(
    pathways = gene_sets,
    stats = salpha,
    scoreType = "pos",
    minSize = min_size,
    maxSize = max_size
  )
}

#' Perform gene set overrepresentation analysis
#' for each bicluster and annotate cells based on
#' the best match.
#'
#' @description
#' per_cluster_goa loads the required gene set, formats it and performs
#' gene overrepresentation analysis for each bicluster in the cabic
#' object.
#'
#' @param cabic A biclustering object of class "caclust"
#'  as obtained from `caclust`.
#' @inheritParams perform_goa
#' @inheritParams load_gene_set
#' @inheritParams format_gene_sets
#'
#' @return
#' A list contain the goa results for each cluster.
#'
#' @export
per_cluster_gsea <- function(
  cadir,
  caobj,
  org,
  set = "CellMarker",
  min_size = 10,
  max_size = 500,
  filter_literature = TRUE
) {
  stopifnot(is(cadir, "cadir"))

  # Ensure that we deal only with clusters consisting of cells and genes.
  suppressWarnings({
    cadir <- rm_monoclusters(cadir)
  })

  # Load gene sets
  gs <- CAbiNet:::load_gene_set(set = set, org = org)
  gene_sets <- format_gene_sets(
    gs,
    min_size = min_size,
    max_size = max_size,
    filter_literature = filter_literature
  )

  # if (isTRUE(filter_literature)) {
  #   gene_sets <- gene_sets[!grepl("et_al\\.", names(gene_sets))]
  # }

  gc <- gene_clusters(cadir)
  gc_names <- sort(unique(gc))

  cadir_rnk <- cadir
  # cadir_rnk@parameters$alpha_genes <- deg2rad(90)
  # Rank all genes
  cadir_rnk <- rank_genes(
    cadir = cadir_rnk,
    caobj = caobj,
    cluster_only = FALSE
  )

  # Perform goa for each cluster
  gsea_res <- list()

  for (c in seq_len(length(gc_names))) {
    cls <- as.character(gc_names[c])
    salpha <- cadir_rnk@gene_ranks[[cls]]$Score
    names(salpha) <- cadir_rnk@gene_ranks[[cls]]$Rowname

    gsea <- perform_gsea(
      salpha = salpha,
      gene_sets = gene_sets,
      min_size = min_size,
      max_size = max_size,
    )
    colnames(gsea)[colnames(gsea) == "pathway"] <- "gene_set"
    gsea <- gsea[order(gsea$padj), ]

    gsea_res[[cls]] <- gsea
  }

  return(gsea_res)
}

#' Perform gene set overrepresentation analysis
#' for each bicluster and annotate cells based on
#' the best match.
#'
#' @description
#' per_cluster_goa loads the required gene set, formats it and performs
#' gene overrepresentation analysis for each bicluster in the cabic
#' object.
#'
#' @param cabic A biclustering object of class "caclust"
#'  as obtained from `caclust`.
#' @inheritParams perform_goa
#' @inheritParams load_gene_set
#'
#' @return
#' A list contain the goa results for each cluster.
#'
#' @export
per_cluster_goa <- function(
  cabic,
  universe,
  org,
  set = "CellMarker",
  min_size = 10,
  max_size = 500,
  filter_literature = FALSE,
  verbose = TRUE
) {
  stopifnot(is(cabic, "caclust"))

  # Ensure that we deal only with clusters consisting of cells and genes.
  suppressWarnings({
    cabic <- rm_monoclusters(cabic)
  })

  # Load gene sets
  gs <- load_gene_set(set = set, org = org)
  gene_sets <- format_gene_sets(
    gs,
    min_size = min_size,
    max_size = max_size,
    filter_literature = filter_literature
  )

  gc <- gene_clusters(cabic)
  gc_names <- sort(unique(gc))

  gc_list <- lapply(X = gc_names, FUN = function(x) {
    names(gc[which(gc == gc_names[x])])
  })

  names(gc_list) <- gc_names

  # Perform goa for each cluster
  goa_res <- list()

  for (c in seq_len(length(gc_list))) {
    gene <- gc_list[[c]]
    clst_name <- names(gc_list)[c]

    goa <- perform_goa(
      gois = gene,
      gene_sets = gene_sets,
      universe = universe,
      min_size = min_size,
      max_size = max_size,
      verbose = verbose
    )

    goa_res[[clst_name]] <- goa
  }

  return(goa_res)
}

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
  goa_res$log_padj <- log10(goa_res$padj)

  cost_mat <- stats::reshape(
    data = goa_res[, c("cluster", "gene_set", "log_padj")],
    direction = "wide",
    idvar = "cluster",
    timevar = "gene_set",
    new.row.names = seq_len(length(unique(goa_res$cluster)))
  )

  cost_mat[is.na(cost_mat)] <- 0
  colnames(cost_mat) <- gsub("log_padj.", "", colnames(cost_mat), fixed = TRUE)

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
    log_padj = -cost_mat[assignments]
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

    # TODO: Should I use NES for GSEA celltype assignment?

    # Solve assignment problem with the hungarian algorithm.
    if (isFALSE(logpval)) {
      cluster_anno <- assign_cts(goa_res)
    } else if (isTRUE(logpval)) {
      cluster_anno <- assign_cts_logpval(goa_res)
    } else {
      rlang::abort("logpval can only be TRUE or FALSE!")
    }

    # Rename clusters based on GSE.
    for (c in seq_len(length(allcs))) {
      noanno <- allcs[c]

      if (!allcs[c] %in% cluster_anno$cluster) {
        ct <- noanno
      } else {
        ct <- cluster_anno[cluster_anno$cluster == allcs[c], "cell_type"]
        log_padj <- cluster_anno[cluster_anno$cluster == allcs[c], "log_padj"]
        if (log_padj < (-log10(alpha))) ct <- noanno
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
    ...,
    caobj = NULL,
    method = "goa", # TODO: add parameter to all function calls elsewhere.
    filter_literature = TRUE # TODO: add parameter to all function calls elsewhere.
  ) {
    stopifnot(methods::is(obj, "caclust"))

    if (method == "goa") {
      enr_res <- CAbiNet::per_cluster_goa(
        cabic = obj,
        universe = universe,
        set = set,
        org = org,
        min_size = min_size,
        max_size = max_size
      )
    } else if (method == "gsea") {
      if (is.null(caobj)) {
        rlang::abort("Please provide a cacomp object for parameter 'caobj'.")
      }
      enr_res <- per_cluster_gsea(
        cadir = obj,
        caobj = caobj,
        set = set,
        org = org,
        min_size = min_size,
        max_size = max_size,
        filter_literature = filter_literature
      )
    } else {
      rlang::abort("Please pick a valid method: Can be either 'goa' or 'gsea'.")
    }

    cabic <- annotate_by_goa(
      obj = obj,
      goa_res = enr_res,
      alpha = alpha
    )

    return(cabic)
  }
)

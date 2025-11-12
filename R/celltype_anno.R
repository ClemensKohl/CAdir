#' @include classes.R
#' @importFrom CAbiNet annotate_biclustering
#' @importFrom CAbiNet annotate_by_goa
NULL

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
  gene_sets <- CAbiNet:::format_gene_sets(gs)

  if (isTRUE(filter_literature)) {
    gene_sets <- gene_sets[!grepl("et_al\\.", names(gene_sets))]
  }

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
    cluster_anno <- assign_cts(goa_res)

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

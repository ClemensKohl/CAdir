#' @include classes.R

NULL
# TODO: Fix documentation

#' Load the required gene set.
#' @description
#' Loads the speciefied gene set and subsets to the required organism.
#' @param set Name of the gene set. Currently only supports "CellMarker"
#' @param org Short name of the organism. "mm" for mouse, "hs" for human.
#' @returns
#' data frame with columns "cell_type" and "gene".
#' @export
load_ct_gene_set <- function(set = "CellMarker", org) {
  stopifnot(org %in% c("mm", "hs"))
  if (set == "CellMarker") {
    gs <- CAdir::cellmarker_v2

    if (org == "mm") {
      gs <- gs[gs$species == "Mouse", ]
    }
    if (org == "hs") {
      gs <- gs[gs$species == "Human", ]
    }

    gs <- gs[, c("cell_name", "marker")]
  } else {
    stop("Other gene sets besides CellMarker are not implemented yet.")
  }

  colnames(gs) <- c("cell_type", "gene")

  return(gs)
}

#' Changes long format gene set data frame into a named list of gene sets.
#'
#' @description
#' Takes in a data frame with 2 columns (gene set name and gene name)
#' and accumulates all genes in the same gene set into a list element.
#'
#' @param gene_sets data frame with first column gene set names and the
#'  second column the gene symbol. Long format expected.
#' @param filter_literature Filter out gene sets
#' named after publications.
#'
#' @returns
#' A named list with gene sets as names and genes as vectors.
make_gene_set_list <- function(gene_sets, filter_literature = FALSE) {
  # Input gene_sets assumed to be a long format data frame
  # Bring data frame into shape that you can use for hypergeom test

  if (is(gene_sets, "tbl")) {
    gene_sets <- as.data.frame(gene_sets)
  }

  gene_list <- list()

  # unique gene set names
  gs <- sort(unique(drop(gene_sets[, 1])))

  for (i in seq_len(length(gs))) {
    sel <- which(gene_sets[, 1] == gs[i])
    genes <- drop(gene_sets[sel, 2])
    gene_list[[gs[i]]] <- genes
  }

  names(gene_list) <- gsub(" ", "_", names(gene_list))

  if (isTRUE(filter_literature)) {
    gene_list <- gene_list[!grepl("et_al\\.", names(gene_list))]
  }
  return(gene_list)
}

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
filter_gene_set_list <- function(
  gene_sets,
  min_size = 10,
  max_size = 500
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
  return(gene_sets[idx])
}

#' Adapted from `DOSE:::enricher_interal`.
#' Performs Gene Overrepresentation Analysis.
#'
#' @description
#' compute_goa takes a number of genes of interest (gois)
#' and a list of gene sets and performs gene overrepresentation
#' analysis.
#'
#' @param gois Genes of interest. Usually the co-clustered genes.
#' @param universe All genes in data set.
#' @param gene_sets Named list of gene sets and their genes.
#' @param verbose Toggles verbosity of warnings.
#' @inheritParams filter_gene_set_list
#'
#' @details
#' compute_goa performs a hypergeometric test on the gois and gene sets.
#' Gene sets are trimmed to genes that are present in universe
#'  (usually all the genes available in the data set) and trimmed to set
#'  min and max size.
#'
#' @returns
#' A data frame with the results of the go-analysis.
#'
#' @export
compute_goa <- function(
  gois,
  gene_sets,
  universe,
  min_size,
  max_size,
  verbose = TRUE
) {
  # subset gene sets to genes in universe
  gene_sets <- lapply(gene_sets, intersect, universe)

  # Subset gois to genes in gene sets
  all_gs <- unique(unlist(gene_sets))
  gois <- gois[gois %in% all_gs]

  # Filter out very small and very large gene sets
  # We do this after subsetting the gois.
  gene_sets <- filter_gene_set_list(
    gene_sets = gene_sets,
    min_size = min_size, # 10
    max_size = max_size
  ) # 500
  # number of clustered genes in each gene set
  gois_in_set <- sapply(gene_sets, intersect, gois)

  # Remove gene sets with 0 gois in them
  gois_in_set <- filter_gene_set_list(
    gene_sets = gois_in_set,
    min_size = 1,
    max_size = Inf
  )
  if (length(gois_in_set) == 0) {
    if (isTRUE(verbose)) {
      warning(
        "No genes of interest are found in any gene set!",
        " Returning empty results data.frame."
      )
    }
    enrich_res <- data.frame(
      gene_set = "No_Cell_Type_Found",
      pval = NA,
      padj = NA,
      GeneRatio = NA,
      BgRatio = NA,
      ngois_in_set = 0,
      ngenes_in_set = length(gene_sets),
      ngois = length(gois),
      ngenes_in_sets = length(all_gs)
    )
    return(enrich_res)
  } else {
    # subset gene sets to those with gois in them
    # make sure the two sets are the same order.
    gs_names <- sort(unique(names(gois_in_set)))
    gene_sets <- gene_sets[gs_names]
    gois_in_set <- gois_in_set[gs_names]

    # Build parameter data frame
    ngois <- length(gois) # clustered genes
    group1 <- lengths(gene_sets) # the size of gene sets
    group2 <- length(all_gs) # total genes in gene sets
    overlap <- lengths(gois_in_set) # nr gois in gene_set

    phyper_df <- data.frame(
      gois_in_set = overlap - 1, # white balls drawn / gois in gs
      genes_in_set = group1, # total white balls / genes in gs
      genes_universe = group2 - group1, # total black balls / all genes in gene set - gs
      ngois = ngois # balls drawn / number gois
    )
    rownames(phyper_df) <- names(gene_sets)

    # Hypergeometric test
    pvalues <- apply(phyper_df, 1, function(n) {
      stats::phyper(n[1], n[2], n[3], n[4], lower.tail = FALSE)
    })

    # adjust p-values
    p_adj <- stats::p.adjust(pvalues, method = "BH")

    ## gene ratio and background ratio
    gene_ratio <- apply(data.frame(a = overlap, b = ngois), 1, function(x) {
      paste(x[1], "/", x[2], sep = "", collapse = "")
    })

    bg_ratio <- apply(data.frame(a = group1, b = group2), 1, function(x) {
      paste(x[1], "/", x[2], sep = "", collapse = "")
    })

    # return results either as data frame or list
    enrich_res <- data.frame(
      gene_set = names(gene_sets),
      pval = pvalues,
      padj = p_adj,
      GeneRatio = gene_ratio,
      BgRatio = bg_ratio,
      ngois_in_set = overlap,
      ngenes_in_set = group1,
      ngois = ngois,
      ngenes_in_sets = group2
    )

    rownames(enrich_res) <- NULL

    ord <- order(enrich_res$padj)

    return(enrich_res[ord, ])
  }
}

#' Perform GSEA using fgsea.
#' @inheritParams compute_goa
#' @param salpha Named vector containing the Salpha gene scores.
#' @returns
#' Data frame of enriched gene sets.
compute_gsea <- function(
  salpha,
  gene_sets,
  min_size = 15,
  max_size = 500
) {
  gsea <- fgsea::fgsea(
    pathways = gene_sets,
    stats = salpha,
    scoreType = "pos",
    minSize = min_size,
    maxSize = max_size
  )
  return(gsea)
}

#' Perform gene set overrepresentation analysis
#' for each bicluster and annotate cells based on
#' the best match.
#'
#' @description
#' gsea_per_cluster loads the required gene set, formats it and performs
#' gene overrepresentation analysis for each bicluster in the cadir
#' object.
#'
#' @param cadir A biclustering object of class "cadir".
#' @param caobj A CA object as obtained from APL.
#' @inheritParams compute_gsea
#' @inheritParams load_ct_gene_set
#' @inheritParams make_gene_set_list
#'
#' @return
#' A list contain the goa results for each cluster.
#'
#' @export
gsea_per_cluster <- function(
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
  gs <- load_ct_gene_set(set = set, org = org)
  gene_sets <- make_gene_set_list(
    gs,
    filter_literature = filter_literature
  )

  gc <- cadir@gene_clusters
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

    gsea <- compute_gsea(
      salpha = salpha,
      gene_sets = gene_sets,
      min_size = min_size,
      max_size = max_size
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
#' goa_per_cluster loads the required gene set, formats it and performs
#' gene overrepresentation analysis for each bicluster in the cadir
#' object.
#'
#' @param cadir A biclustering object of class "cadir".
#' @inheritParams compute_goa
#' @inheritParams load_ct_gene_set
#' @inheritParams make_gene_set_list
#'
#' @return
#' A list contain the goa results for each cluster.
#'
#' @export
goa_per_cluster <- function(
  cadir,
  universe,
  org,
  set = "CellMarker",
  min_size = 10,
  max_size = 500,
  filter_literature = FALSE,
  verbose = TRUE
) {
  stopifnot(is(cadir, "cadir"))

  # Ensure that we deal only with clusters consisting of cells and genes.
  suppressWarnings({
    cadir <- rm_monoclusters(cadir)
  })

  # Load gene sets
  gs <- load_ct_gene_set(set = set, org = org)
  gene_sets <- make_gene_set_list(
    gs,
    filter_literature = filter_literature
  )

  gc <- cadir@gene_clusters
  gc_names <- sort(unique(gc))

  gc_list <- lapply(X = gc_names, FUN = function(x) {
    names(gc[which(gc == gc_names[x])])
  })

  names(gc_list) <- gc_names

  # Perform goa for each cluster
  gse_res <- list()

  for (c in seq_len(length(gc_list))) {
    gene <- gc_list[[c]]
    clst_name <- names(gc_list)[c]

    goa <- compute_goa(
      gois = gene,
      gene_sets = gene_sets,
      universe = universe,
      min_size = min_size,
      max_size = max_size,
      verbose = verbose
    )

    gse_res[[clst_name]] <- goa
  }

  return(gse_res)
}

#' Assign cell types to clusters using the Hungarian algorithm.
#' @description
#' Uses the hungarian algorithm (assignment problem)
#' to assign a cell type from the gene set overrepresentation
#' analysis to one (and only one) cluster.
#' The cost matrix is based on the adjusted p-value.
#' @param gse_res List of gene set enrichment results for each bicluster.
#' @param cost Column to use for the cost matrix. Can be
#' "pval", "logpval" or "NES".
#' @returns
#' A data frame with the assigned cell types and adjusted p-values.
#' @export
run_hungarian <- function(gse_res, cost = "pval") {
  # Solve assignment problem with the hungarian algorithm.
  # Build cost matrix.
  gse_res <- dplyr::bind_rows(gse_res, .id = "cluster")

  if (cost == "pval") {
    gse_res$cost <- gse_res[, "padj"]
    na_cost <- 1
    result_factor <- 1
  } else if (cost == "logpval") {
    gse_res$cost <- log10(gse_res$padj)
    na_cost <- 0
    result_factor <- -1
  } else if (cost == "NES") {
    gse_res$cost <- -gse_res$NES
    na_cost <- 0
    result_factor <- -1
  } else {
    rlang::abort("Invalid cost.")
  }

  cost_mat <- stats::reshape(
    data = gse_res[, c("cluster", "gene_set", "cost")],
    direction = "wide",
    idvar = "cluster",
    timevar = "gene_set",
    new.row.names = seq_len(length(unique(gse_res$cluster)))
  )

  cost_mat[is.na(cost_mat)] <- na_cost
  colnames(cost_mat) <- gsub("cost.", "", colnames(cost_mat), fixed = TRUE)

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
    cost = result_factor * cost_mat[assignments]
  )

  return(cluster_anno)
}


#' Annotate CAdir results by gene set enrichment
#' @description
#' `annotate_by_gse` takes a biclustering results
#' and annotates it with the gene overrepresentation analysis results (goa).
#'
#' @param gse_res List of goa results for each bicluster.
#' @param obj `cadir` object with biclustering results. Alternatively could be
#' a caclust object.
#' @param cost value that should be
#' used for constructing the assignment cost matrix.
#' Can bei either "pval" (adjusted p-value),
#' "logpval" (adjusted log10(p-value)) or "NES" (normalized enrichment score,
#' only applicable for GSEA)
#' @param alpha p-value cutoff for annotating clusters when using cost "pval" or "logpval"
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
annotate_by_gse <- function(
  obj,
  gse_res,
  alpha = 0.05,
  cost = c("pval", "logpval", "NES")
) {
  stopifnot(methods::is(obj, "cadir"))
  if (length(cost) > 1) {
    cost <- cost[1]
  }
  if (!cost %in% c("pval", "logpval", "NES")) {
    rlang::abort("Please choose a valid 'cost'!")
  }
  # dictionary
  dict <- obj@dict

  # cell clusters
  ccs <- obj@cell_clusters
  ccs_nm <- names(ccs)
  unccs <- levels(ccs)

  # gene clusters
  gcs <- obj@gene_clusters
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
  nclust <- length(gse_res)

  gse_res <- lapply(gse_res, function(x, nc = nclust) {
    subs <- min(nrow(x), nc)
    x[seq_len(subs), ]
  })

  # Solve assignment problem with the hungarian algorithm.
  cluster_anno <- run_hungarian(gse_res, cost = cost)

  # Rename clusters based on GSE.
  for (c in seq_len(length(allcs))) {
    noanno <- allcs[c]

    if (!allcs[c] %in% cluster_anno$cluster) {
      ct <- noanno
    } else {
      ct <- cluster_anno[cluster_anno$cluster == allcs[c], "cell_type"]
      cost_val <- cluster_anno[cluster_anno$cluster == allcs[c], "cost"]

      if (cost == "pval") {
        if (cost_val > alpha) ct <- noanno
      } else if (cost == "logpval") {
        if (cost_val < (-log10(alpha))) ct <- noanno
      } else if (cost == "NES") {
        # if (cost_val < (-log10(alpha))) ct <- noanno
      } else {
        rlang::abort()
      }
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

#' Perform gene overrepresentation analysis and annotate biclusters.
#'
#' @description
#' Wrapper function for `goa_per_cluster` and `annotate_by_goa`.
#'
#' @param obj Either a `cadir`.
#' @inheritParams goa_per_cluster
#' @inheritParams gsea_per_cluster
#' @inheritParams annotate_by_gse
#' @param method Either "goa" for gene overrepresentation analysis or "gsea"
#' for gene set enrichment analysis.
#'
#' @details
#' `annotate` performs per cluster GOA with a hypergeometric
#'  and annotates the biclustering results from CAdir.
#'
#' @returns
#' An object of type `cadir` with annotated biclusters.
#' @export
setGeneric(
  "annotate",
  function(
    obj,
    universe,
    org,
    caobj = NULL,
    set = "CellMarker",
    alpha = 0.05,
    min_size = 10,
    max_size = 500,
    verbose = TRUE,
    method = "goa", # TODO: add parameter to all function calls elsewhere.
    filter_literature = FALSE, # TODO: add parameter to all function calls elsewhere.
    cost = c("pval", "logpval", "NES")
  ) {
    standardGeneric("annotate")
  }
)

#' Annotate the biclustering
#' @inheritParams annotate
#' @inheritParams annotate_by_gse
#' @param obj A cadir object.
#' @rdname annotate
#' @export
setMethod(
  f = "annotate",
  signature = (obj <- "cadir"),
  function(
    obj,
    universe,
    org,
    caobj = NULL,
    set = "CellMarker",
    alpha = 0.05,
    min_size = 10,
    max_size = 500,
    verbose = TRUE,
    method = "goa",
    filter_literature = FALSE,
    cost = c("pval", "logpval", "NES")
  ) {
    stopifnot(methods::is(obj, "cadir"))

    if (length(cost) > 1) {
      cost <- cost[1]
    }

    if (!cost %in% c("pval", "logpval", "NES")) {
      rlang::abort("Please choose a valid 'cost'!")
    }

    if (method != "gsea" && cost == "NES") {
      rlang::abort("NES cost can only be used with GSEA.")
    }

    if (method == "goa") {
      enr_res <- goa_per_cluster(
        cadir = obj,
        universe = universe,
        set = set,
        org = org,
        min_size = min_size,
        max_size = max_size,
        filter_literature = filter_literature
      )
    } else if (method == "gsea") {
      if (is.null(caobj)) {
        rlang::abort("Please provide a cacomp object for parameter 'caobj'.")
      }
      enr_res <- gsea_per_cluster(
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

    cadir <- annotate_by_gse(
      obj = obj,
      gse_res = enr_res,
      alpha = alpha,
      cost = cost
    )

    return(cadir)
  }
)

# NOTE: Verbosity of the messages can be controlled with rlang.
# To turn all messages off:
# options(rlib_message_verbosity = "quiet")
# # To turn them back on:
# options(rlib_message_verbosity = "default")

#' Split a cluster into 2 sub-clusters.
#' @inheritParams dirclust
#' @param idx The cluster index to split.
#' @param cadir A `cadir` object.
#' @returns
#' A `cadir` object that contains only the two new clusters.
sub_cluster <- function(cadir, points, idx) {
  stopifnot(idx %in% levels(cadir@cell_clusters))

  sel <- which(cadir@cell_clusters == idx)
  cl <- points[sel, ]
  if (length(cl) == 0) {
    rlang::abort("No points in cluster")
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
  min_cells = 5,
  make_plots = FALSE
) {
  cls <- levels(cadir@cell_clusters)

  for (i in cls) {
    stopifnot(i %in% levels(cadir@cell_clusters))

    if (sum(cadir@cell_clusters == i) < 2) {
      next
    }

    sres <- sub_cluster(
      cadir = cadir,
      points = caobj@prin_coords_cols,
      idx = i
    )

    if (length(unique(sres@cell_clusters)) < 2) {
      next
    }
    stopifnot(length(unique(sres@cell_clusters)) == 2)

    elems <- as.numeric(table(sres@cell_clusters))
    if (any(elems < min_cells)) {
      next
    }

    to_split <- decide_split(sres@directions, cutoff = cutoff)

    if (isTRUE(to_split)) {
      rlang::inform(paste0("\tSplitting cluster ", i))

      if (isTRUE(make_plots)) {
        p <- cluster_apl(
          caobj = caobj,
          cadir = sres,
          direction = cadir@directions[cadir@dict[[f2c(i)]], ],
          group = which(cadir@cell_clusters == i),
          cluster = NULL,
          highlight_cluster = FALSE,
          show_genes = FALSE,
          show_lines = TRUE
        )

        rep <- paste0("rep_", cadir@log$last_rep)
        nms_rep <- names(cadir@plots$splits[[rep]])
        nm <- f2c(i)

        if (nm %in% names(nms_rep)) {
          cnt <- sum(grepl(nm, nms_rep))
          nm <- paste0(nm, cnt + 1)
        }

        cadir@plots$splits[[rep]][[nm]] <- p
      }

      # Update the clustering results and add the new directions.
      new_cl <- sort(unique(sres@cell_clusters))

      lvls <- levels(cadir@cell_clusters)
      extra_nr <- length(lvls) + 1
      extra <- cl2nm(extra_nr)

      if (extra %in% lvls) {
        extra_nr <- max(get_std_num(lvls[is_std_name(lvls)]))
        extra <- cl2nm(extra_nr + 1)
      }

      # Just making sure that we arent overwriting an existing cluster.
      # Should not be possible though.
      stopifnot(!extra %in% lvls)

      sel <- which(cadir@cell_clusters == i)
      cl1 <- sres@cell_clusters == new_cl[1]
      cl2 <- sres@cell_clusters == new_cl[2]

      if (!length(cl1) > 1 && length(cl2) > 1) {
        next
      }

      # Update the changed clusters in original cadir:
      # 1) update the factor levels
      new_lvls <- c(lvls, extra)
      levels(sres@cell_clusters) <- new_lvls
      levels(cadir@cell_clusters) <- new_lvls

      # 2) Update the cluster assignments
      sres@cell_clusters[cl1] <- i
      sres@cell_clusters[cl2] <- extra
      cadir@cell_clusters[sel] <- sres@cell_clusters

      # 3) Update the directions
      rownames(sres@directions) <- c(i, extra)
      cadir@directions[i, ] <- sres@directions[i, ]

      rn <- rownames(cadir@directions)
      cadir@directions <- rbind(
        cadir@directions,
        sres@directions[extra, ]
      )

      rownames(cadir@directions) <- c(rn, extra)

      # 4) Add new cluster to dict.
      cadir@dict[[extra]] <- extra_nr

      # 5) distances are wrong. reset.
      cadir@distances <- matrix(0, 0, 0)

      # 6) Make sure we have coherent naming.
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
merge_clusters <- function(caobj, cadir, cutoff, make_plots = FALSE) {
  samples <- caobj@prin_coords_cols
  clusters <- cadir@cell_clusters
  directions <- cadir@directions
  dir_nms <- rownames(directions)

  asim_cutoff <- 1 - cutoff / pi
  sim <- get_ang_sim(directions, directions)
  # Set the lower diagonal to 0
  sim[lower.tri(sim, diag = TRUE)] <- 0
  candidates <- apply(
    X = sim,
    MARGIN = 1,
    FUN = function(x) which(x >= asim_cutoff),
    simplify = FALSE
  )

  sel <- which(lengths(candidates) > 0)

  for (s in sel) {
    merge_parent <- search_dict(cadir@dict, s)

    cds <- candidates[[s]]
    cds_nms <- dir_nms[cds]
    merger_idxs <- c(s, cds)

    # get the names of the mergers
    mergers <- search_dict(cadir@dict, merger_idxs)

    rlang::inform(paste0("\tMerging ", merge_parent, " with ", cds_nms))

    # clusters that have to be changed.
    cls <- which(clusters %in% mergers)

    if (length(cls) == 0) {
      next
    }
    if (length(unique(clusters[cls])) == 1) {
      next
    }

    new_dir <- drop(total_least_squares(samples[cls, , drop = FALSE]))

    if (isTRUE(make_plots)) {
      rep <- paste0("rep_", cadir@log$last_rep)

      sres <- methods::new(
        "cadir",
        cell_clusters = droplevels(clusters[cls]),
        directions = directions[mergers, ]
      )

      names(sres@cell_clusters) <- names(clusters[cls])

      p <- cluster_apl(
        caobj = caobj,
        cadir = sres,
        direction = new_dir,
        group = cls,
        cluster = merge_parent,
        highlight_cluster = FALSE,
        show_cells = TRUE,
        show_genes = FALSE,
        show_lines = TRUE
      )

      nms_rep <- names(cadir@plots$merges[[rep]])
      nm <- paste(
        "cluster",
        merger_idxs,
        collapse = "_",
        sep = ""
      )

      if (nm %in% names(nms_rep)) {
        cnt <- sum(grepl(nm, nms_rep))
        nm <- paste0(nm, cnt + 1)
      }

      cadir@plots$merges[[rep]][[nm]] <- p
    }

    # Update orig. cadir with new cluster:
    # 1) Update cluster assignments
    clusters[cls] <- merge_parent
    cadir@cell_clusters <- clusters

    # 2) Update directions
    directions[s, ] <- new_dir
    directions <- directions[-cds, , drop = FALSE]
    cadir@directions <- directions

    # 3) distances arent correct anymore. So set empty.
    cadir@distances <- matrix(0, 0, 0) # dists not true anymore.

    # 4) Ensure dict is up to date
    cadir@dict[[merge_parent]] <- s
    cadir@dict[!names(cadir@dict) %in% cds_nms]

    # 5) Ensure coherent naming.
    cadir <- rename_clusters(cadir = cadir)

    out <- merge_clusters(
      caobj = caobj,
      cadir = cadir,
      cutoff = cutoff,
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
dirclust_splitmerge <- function(
  caobj,
  k,
  cutoff = NULL,
  qcutoff = 0.8,
  method = "random",
  counts = NULL,
  apl_quant = 0.99,
  min_cells = 5,
  epochs = NULL,
  reps = NULL,
  apl_cutoff_reps = 100,
  make_plots = FALSE,
  convergence_thr = 0.001,
  max_iter = 50,
  init = "kmeanspp"
) {
  fun_args <- match.call()

  # Convert cutoff to radians
  if (!is.null(cutoff)) {
    cutoff <- deg2rad(cutoff)
  }

  #########
  # Setup #
  #########

  log <- list()
  log$clusters <- as.data.frame(matrix(
    0,
    ncol = 1,
    nrow = nrow(caobj@prin_coords_cols)
  ))

  colnames(log$clusters) <- "root"

  log$directions <- as.data.frame(matrix(
    0,
    ncol = ncol(caobj@prin_coords_cols),
    nrow = 1
  ))

  colnames(log$directions) <- colnames(caobj@prin_coords_cols)
  log$directions <- cbind(iter = "root", dirname = "root", log$directions)

  ######################
  # Initial clustering #
  ######################

  out <- dirclust(
    points = caobj@prin_coords_cols,
    k = k,
    epochs = epochs,
    init = init,
    lines = NULL,
    log = FALSE,
    convergence_thr = convergence_thr,
    max_iter = max_iter
  )

  out <- rename_clusters(out)

  log <- log_iter(log = log, cadir = out, name = "iter_0")

  ##########
  # Cutoff #
  ##########

  # If the cutoff is set to NULL we computer the S-alpha cutoff.
  if (is.null(cutoff)) {
    if (method == "permutation" && is.null(counts)) {
      rlang::warn(paste0(
        "No count matrix for permutation supplied. ",
        "Switching to random directions method."
      ))
      counts <- NULL
      method <- "random"
    }

    # HACK: Better way than random sampling?
    # We need a group for the permutation method.
    # We just pick a random cluster.
    grp_idx <- which(
      out@cell_clusters == base::sample(unique(out@cell_clusters), 1)
    )

    cutoff <- get_apl_cutoff(
      caobj = caobj,
      counts = counts,
      method = method,
      group = grp_idx,
      quant = apl_quant,
      apl_cutoff_reps = apl_cutoff_reps,
      axis = "cols"
    )

    out@parameters$sa_cutoff <- cutoff

    rlang::inform(paste0(
      "\nInferred cutoff angle: ",
      round(rad2deg(cutoff), 2),
      "\n"
    ))
  }

  use_conv <- is.null(reps)
  if (is.null(reps)) {
    reps <- max_iter
  }
  i <- 0
  has_conv <- FALSE

  while (!(isTRUE(has_conv) || i >= reps)) {
    i <- i + 1
    rlang::inform(paste0("Iteration ", i))

    iter_nm <- paste0("iter_", i)
    out@log$last_rep <- i

    ##################
    # Split clusters #
    ##################

    out <- split_clusters(
      cadir = out,
      caobj = caobj,
      cutoff = cutoff,
      min_cells = min_cells,
      make_plots = make_plots
    )

    log <- log_iter(log = log, cadir = out, name = paste0("split_", iter_nm))

    out <- dirclust(
      points = caobj@prin_coords_cols,
      k = ncol(out@directions),
      epochs = epochs,
      log = FALSE,
      cadir = out,
      convergence_thr = convergence_thr,
      max_iter = 10
    )

    out <- rename_clusters(out)

    log <- log_iter(log = log, cadir = out, name = paste0("interS_", iter_nm))

    ##################
    # Merge clusters #
    ##################

    out <- merge_clusters(
      caobj = caobj,
      cadir = out,
      cutoff = cutoff,
      make_plots = make_plots
    )

    log <- log_iter(log = log, cadir = out, name = paste0("merge_", iter_nm))

    out <- dirclust(
      points = caobj@prin_coords_cols,
      k = ncol(out@directions),
      epochs = epochs,
      log = FALSE,
      cadir = out,
      convergence_thr = convergence_thr,
      max_iter = 10
    )
    out <- rename_clusters(out)

    log <- log_iter(log = log, cadir = out, name = paste0("interM_", iter_nm))

    if (i > 1 && use_conv) {
      has_conv <- check_convergence(
        now = out@directions,
        prev = last_iter,
        convergence_thr = convergence_thr
      )
      if (isTRUE(has_conv)) break
    }

    last_iter <- out@directions
  }

  ################
  # Assign genes #
  ################

  out@gene_clusters <- assign_genes(
    caobj = caobj,
    cadir = out,
    qcutoff = qcutoff,
    coords = "prin"
  )

  out <- rename_clusters(out)

  log <- log_iter(log = log, cadir = out, name = paste0("end", iter_nm))

  ################
  # add log info #
  ################

  to_keep <- (!names(out@log) %in% names(log))
  out@log <- c(out@log[to_keep], log)
  out@parameters$qcutoff <- qcutoff
  out@parameters$call <- fun_args

  return(out)
}

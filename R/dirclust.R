#' Initialize directions by kmeans++ method
#' @param points Row-wise matrix of points to be clustered.
#' @param k Number of clusters.
#' @returns
#' Vector of indices of the initial directions.
kmeanspp_init <- function(points, k) {
  pnorm <- row_norm(points)
  npoints <- nrow(points)

  center_ids <- integer(k)
  center_ids[1] <- sample.int(npoints, 1)

  for (ii in 2:k) {
    lines <- points[center_ids, ] / pnorm[center_ids]

    ldist <- dist_to_line(
      points = points,
      lines = lines,
      pnorm = pnorm,
      pos_only = FALSE
    )

    probs <- apply(ldist, 1, min)
    probs[center_ids] <- 0

    center_ids[ii] <- sample.int(npoints, 1, prob = probs)
  }

  return(center_ids)
}

#' Calculate the distance of points to the directions
#' @param points Row-wise matrix of points.
#' @param lines Row-wise matrix of lines/directions.
#' @param pnorm Vector of the norm of the points.
#' @param pos_only If TRUE only considers points with a
#' positive projection onto the line.
#' @returns
#' Matrix of distances of points to the directions.
dist_to_line <- function(points, lines, pnorm, pos_only = TRUE) {
  if (is.matrix(lines)) {
    lines <- t(lines)
  }
  proj <- (points %*% lines)
  dist <- pnorm^2 - proj^2
  dist[dist < 0] <- 0
  dist <- sqrt(dist)

  # Ensure that we only look at positive directions
  if (pos_only) {
    dist[proj < 0] <- Inf
  }

  return(dist)
}

#' Perform total least squares regression on points.
#' @param points Row-wise matrix of points.
#' @returns
#' Vector of the regression line (the first PC).
total_least_squares <- function(points) {
  if (nrow(points) == 1) {
    reg_line <- points / row_norm(points)
  } else {
    suppressWarnings({
      reg_line <- irlba::irlba(points, nv = 1, right_only = TRUE)$v
    })

    is_flipped <- sign_flip(points = points, line = reg_line)

    # Ensure that line is pointing towards the majority of points.
    if (isTRUE(is_flipped)) {
      reg_line <- reg_line * (-1)
    }
  }

  return(reg_line)
}

#' Update directions through Total Least Squares Regression
#' @param points Row-wise matrix of points.
#' @param clusters Vector of cluster assignments for the points.
#' @param lines Row-wise matrix of lines.
#' @param k Number of clusters.
#' @returns
#' Matrix of updated directions.
update_line <- function(points, clusters, lines, k) {
  # SVD gives total least squares results
  k <- sort(unique(clusters))
  for (c in seq_len(length(k))) {
    sel <- which(clusters == k[c])

    if (length(sel) == 0) {
      next
    } else {
      lines[c, ] <- total_least_squares(points[sel, , drop = FALSE])
    }
  }
  return(lines)
}

#' Basic clustering by directions in CA space
#' @inheritParams check_convergence
#' @param points Matrix of points to be clustered.
#'  Points are expected to be rows.
#' @param k Number of clusters.
#' @param epochs Number of iterations.
#' @param init Initialization method for the lines. Options are 'rand' and
#' 'kmeanspp'.
#' @param lines Optional. Row-wise matrix of lines.
#'  If not NULL, then the provided lines are used for initialization.
#' @param cadir If a cadir object is provided, uses this instead of `lines` to
#' initialize the clustering.
#' @param log Logical. If TRUE, then the function returns a list with
#' the computed distances and directions at each step.
#' @param max_iter Maximum number of iterations if no convergence.
#'
#' @returns
#' An element of class cadir.
#'
#' @export
dirclust <- function(
  points,
  k,
  epochs = NULL,
  init = "kmeanspp",
  lines = NULL,
  log = FALSE,
  cadir = NULL,
  convergence_thr = 0.001,
  max_iter = 50
) {
  if (!is.null(cadir)) {
    stopifnot(methods::is(cadir, "cadir"))
    if (!is.null(lines)) {
      rlang::warn(paste0(
        "You overspecified the init. directions ",
        "with 'lines' and 'cadir'. ",
        "Picking directions from 'cadir'."
      ))
    }
    lines <- cadir@directions
  }

  pnorm <- row_norm(points)

  if (is.null(lines)) {
    if (init == "rand") {
      # initialize line
      idxs <- rand_idx(points, k)
      lines <- points[idxs, ]
      lines <- lines / row_norm(lines)
    } else if (init == "kmeanspp") {
      center_ids <- kmeanspp_init(points = points, k = k)
      lines <- points[center_ids, ]
    } else {
      rlang::abort("Please use init = 'rand' or 'kmeanspp'.")
    }
  } else {
    # ensure that they are unit vectors.
    lines <- lines / row_norm(lines)
  }

  rownames(lines) <- paste0("cluster_", seq_len(nrow(lines)))

  if (isTRUE(log)) {
    dist_log <- vector(mode = "list", length = epochs)

    dir_log <- vector(mode = "list", length = epochs)
    names(dir_log) <- as.character(seq_len(length(dir_log)))
    dir_log[["0"]] <- lines
  }

  use_conv <- is.null(epochs)
  if (is.null(epochs)) {
    epochs <- max_iter
  }
  i <- 0
  has_conv <- FALSE

  while (!(isTRUE(has_conv) || i >= epochs)) {
    i <- i + 1

    # calculate distance to line.
    ldist <- dist_to_line(
      points = points,
      lines = lines,
      pnorm = pnorm,
      pos_only = TRUE
    )

    # find closest line
    clusters <- apply(ldist, 1, which.min)

    # update lines
    lines <- update_line(points, clusters, lines, k)

    if (isTRUE(log)) {
      dir_log[[as.character(i)]] <- lines

      dist_log[[i]] <- vapply(
        seq_len(nrow(ldist)),
        function(x) {
          ldist[x, clusters[x]]
        },
        numeric(1)
      )
    }

    if (i > 1 && use_conv) {
      has_conv <- check_convergence(
        now = lines,
        prev = last_iter,
        convergence_thr = convergence_thr
      )
      if (isTRUE(has_conv)) break
    }

    last_iter <- lines
  }

  if (isTRUE(log)) {
    dist_log <- as.data.frame(do.call(rbind, dist_log))
    dimnames(dist_log) <- list(
      paste0("epoch_", seq_len(nrow(dist_log))),
      rownames(points)
    )
    dist_log$epoch <- seq_len(nrow(dist_log))
    log_list <- list("dist_log" = dist_log, "dir_log" = dir_log)
  } else {
    log_list <- list()
  }

  # Remove lines without clusters
  uni_clust <- sort(unique(clusters))
  if (length(uni_clust) != nrow(lines)) {
    lines <- lines[uni_clust, , drop = FALSE]
    ldist <- ldist[, uni_clust, drop = FALSE]
  }

  cnms <- names(clusters)
  clusters <- paste0("cluster_", clusters)
  clusters <- factor(clusters, levels = sort(unique(clusters)))
  names(clusters) <- cnms

  dict <- as.list(seq_len(nrow(lines)))
  names(dict) <- rownames(lines)

  if (is.null(cadir)) {
    out <- methods::new(
      "cadir",
      cell_clusters = clusters,
      directions = lines,
      distances = ldist,
      parameters = list(
        "k" = k,
        "epochs" = epochs,
        "init" = init,
        "log" = log
      ),
      log = log_list,
      dict = dict
    )
  } else {
    out <- cadir
    out@cell_clusters <- clusters
    out@directions <- lines
    out@distances <- ldist
    to_keep <- (!names(out@log) %in% names(log))
    out@log <- c(out@log[to_keep], out@log)
    out@dict <- dict
  }

  return(out)
}

#' Assigns cells to the closest direction.
#' @param cells Row-wise matrix of cell coordinates.
#' @param directions Row-wise matrix of directions.
#' @returns
#' Vector of same length as nrow(cells) with row-indices of the direction the
#' cell is assigned to.
#' @export
assign_cells <- function(cells, directions) {
  # vector norm
  pnorm <- row_norm(cells)

  # calculate distance to line.
  ldist <- dist_to_line(cells, directions, pnorm, pos_only = TRUE)

  # find closest line
  clusters <- apply(ldist, 1, which.min)
  clusters <- rownames(directions)[clusters]

  return(clusters)
}


#' Determine sign for SVD singular vectors.
#'
#' @param points Points/Cells for which to check the directionality of the line.
#' @param line Line/direction whose directionality needs to be determined.
#'
#' @references
#' Bro, Rasmus, Acar, Evrim, and Kolda, Tamara Gibson.
#' Resolving the sign ambiguity in the singular value decomposition.
#' United States: N. p., 2007. Web. doi:10.2172/920802.
#'
#' @returns
#' TRUE if the sign of the direction needs to be flipped, FALSE otherwise.
sign_flip <- function(points, line) {
  s <- sum(sign(points %*% line) * (points %*% line)**2)
  return(s < 0)
}

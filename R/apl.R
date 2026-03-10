#' Create a model to project points into an Association Plot
#' @param coords The CA coordinates to use (should be standard or principal coordinates).
#' @param direction Normed direction vector of the APL plot (usually same coordinate type as coords).
#' @param group A vector of indices which indicate the points
#' that belong to the cluster. Only needed here to orient the plot.
#' @returns
#' A model that can be used to project new points onto the APL plot.
#' @export
apl_model <- function(
  coords,
  direction,
  group = NULL
) {
  # stopifnot(methods::is(caobj, "cacomp"))

  # coords <- caobj@prin_coords_cols

  avg_group_coords <- direction
  length_vector_group <- sqrt(drop(avg_group_coords %*% avg_group_coords))

  # The line sometimes point into the "wrong" direction.
  # We can determine the cosine, if its negative we flip the x coords.

  cosangle <- 1

  if (!is.null(group)) {
    subgroup <- coords[group, ]

    if (length(group) == 1) {
      group_mean <- subgroup # single sample
    } else {
      group_mean <- colMeans(subgroup) # centroid vector.
    }

    cosangle <- cosine(a = group_mean, b = direction)

    # group_norm <- row_norm(group_mean)
    # gx <- drop(group_mean %*% avg_group_coords) / group_norm
    #
    # if (sign(gx) == -1) {
    #     avg_group_coords <- -avg_group_coords
    # }

    if (cosangle < 0) {
      avg_group_coords <- -avg_group_coords
    }
  }

  model <- function(vec) {
    length_vector <- row_norm(vec)
    cordx <- drop(vec %*% avg_group_coords) / length_vector_group
    cordy <- suppressWarnings(sqrt(length_vector^2 - cordx^2))

    cordx[is.na(cordx)] <- 0
    cordy[is.na(cordy)] <- 0

    return(cbind("x" = cordx, "y" = cordy))
  }

  return(model)
}


#' Random direction association plot coordinates
#'
#' @description
#' Calculates cotangent values for random directions
#'
#' @inheritParams get_apl_cutoff
#' @param dims Number of dimensions to use for the analysis.
#'
#' @returns
#' Numeric vector of cotangent values of length `row_num * apl_cutoff_reps`.
random_direction_cutoff <- function(
  caobj,
  dims = caobj@dims,
  apl_cutoff_reps = 300,
  axis = "cols"
) {
  if (axis == "cols") {
    cols <- caobj@prin_coords_cols
    ca_coords <- caobj@prin_coords_cols
  } else if (axis == "rows") {
    cols <- caobj@std_coords_cols
    ca_coords <- caobj@prin_coords_rows
  }

  row_num <- nrow(ca_coords)

  if (caobj@dims == 1 && !is.empty(caobj@dims)) {
    row_num <- 1
  }

  if (dims < caobj@dims) {
    caobj <- APL::subset_dims(caobj = caobj, dims = dims)
  }

  cols <- caobj@prin_coords_cols
  length_vector_coords <- sqrt(rowSums(ca_coords^2))
  apl_cotan <- numeric(row_num * apl_cutoff_reps)

  for (k in seq(apl_cutoff_reps)) {
    avg_group_coords <- stats::runif(
      n = dims,
      min = 0,
      max = stats::quantile(cols, 0.99)
    )

    length_vector_group <- sqrt(drop(avg_group_coords %*% avg_group_coords))
    coord_x <- drop(ca_coords %*% avg_group_coords) / length_vector_group
    coord_y <- sqrt(pmax(length_vector_coords^2 - coord_x^2, 0))
    cotan <- coord_x / coord_y
    cotan[is.na(cotan)] <- 0

    idx <- ((1:row_num) + ((k - 1) * row_num))
    apl_cotan[idx] <- cotan
  }
  return(apl_cotan)
}


#' Calculates permuted association plot coordinates
#'
#' @description
#' Calculates matrix of apl coordinates when permuting the original data.
#'
#' @inheritParams get_apl_cutoff
#' @param python If TRUE, use python for CA SVD.
#' @param dims The number of dimensions to when performing SVD.
#' Usually can be kept as `caobj@dims`.
#' @returns
#' Numeric vector of cotangent values of length `row_num * apl_cutoff_reps`.
#'
permutation_cutoff <- function(
  caobj,
  counts,
  group = caobj@group,
  dims = caobj@dims,
  apl_cutoff_reps = 10,
  python = TRUE,
  axis = "cols"
) {
  if (axis == "cols") {
    row_num <- nrow(caobj@V)
  } else if (axis == "rows") {
    row_num <- nrow(caobj@U)
  }

  apl_cotan <- numeric(row_num * apl_cutoff_reps)

  if (caobj@dims == 1 && !is.empty(caobj@dims)) {
    row_num <- 1
  }

  margin <- 1
  pc <- 1
  cc <- (axis == "cols")
  cr <- (axis == "rows")

  for (k in seq(apl_cutoff_reps)) {
    mat_perm <- t(apply(counts, margin, FUN = sample))
    colnames(mat_perm) <- colnames(counts)

    suppressWarnings({
      caobjp <- APL::cacomp(
        obj = mat_perm,
        python = python,
        coords = TRUE,
        princ_coords = pc,
        dims = dims,
        top = caobj@top_rows,
        residuals = caobj@params$residuals,
        clip = caobj@params$clip,
        cutoff = caobj@params$cutoff,
        rm_zeros = caobj@params$rm_zeros,
        inertia = FALSE
      )
    })

    caobjp <- APL::apl_coords(
      caobj = caobjp,
      group = group,
      calc_cols = cc,
      calc_rows = cr
    )

    apl_xy <- if (axis == "cols") caobjp@apl_cols else caobjp@apl_rows
    cotan <- apl_xy[, 1] / apl_xy[, 2]
    cotan[is.na(cotan)] <- 0

    idx <- ((seq_len(row_num) + ((k - 1) * row_num)))
    apl_cotan[idx] <- cotan
  }

  return(apl_cotan)
}


#' Calculates the S_alpha cutoff based on random directions or permutations of the data
#'
#' @param caobj A cacomp object.
#' @param method Method to use for computing the cutoff.
#' Either "random" or "permutation".
#' @param group A vector of group indices.
#' @param counts The original count matrix that was used to compute the caobj.
#' @param quant The quantile to use for the cutoff.
#' @param apl_cutoff_reps The number of repetitions to use.
#' Should be between 3-10 for permutation, and >=100 for random.
#' @param axis The axis for which the cutoff should be calculated, either "cols" or "rows".
#'
#' @returns
#' The cutoff angle alpha in radians.
get_apl_cutoff <- function(
  caobj,
  method = "random",
  group = caobj@group,
  counts = NULL,
  quant = 0.99,
  apl_cutoff_reps = 100,
  axis = "cols"
) {
  if (method == "random") {
    if (is.null(apl_cutoff_reps)) {
      apl_cutoff_reps <- 100
    } else if (apl_cutoff_reps < 100) {
      rlang::warn("Number of repetitions should be set >=100.")
      apl_cutoff_reps <- 100
    }

    apl_cotan <- random_direction_cutoff(
      caobj = caobj,
      dims = caobj@dims,
      apl_cutoff_reps = apl_cutoff_reps,
      axis = axis
    )
  } else if (method == "permutation") {
    if (is.null(apl_cutoff_reps)) {
      apl_cutoff_reps <- 5
    } else if (apl_cutoff_reps > 10) {
      rlang::inform("Large number of repetitions might take a long time.")
    }

    apl_cotan <- permutation_cutoff(
      caobj = caobj,
      counts = counts,
      group = group,
      dims = caobj@dims,
      apl_cutoff_reps = apl_cutoff_reps,
      python = TRUE,
      axis = axis
    )
  }

  cutoff_cotan <- stats::quantile(apl_cotan, quant)

  # tan = 1/cotan
  # angle alpha is in radian.
  alpha <- atan(1 / cutoff_cotan)

  return(alpha)
}

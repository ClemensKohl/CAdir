#' Create a model to project points into an Association Plot
#' @param caobj A cacomp object.
#' @param direction Normed direction vector of the APL plot.
#' @param group A vector of indices which indicate the points
#' that belong to the cluster. Only needed here to orient the plot.
#' @returns
#' A model that can be used to project new points onto the APL plot.
#' @export
apl_model <- function(
    caobj,
    direction,
    group = NULL) {
    stopifnot(methods::is(caobj, "cacomp"))

    cent <- caobj@prin_coords_cols

    avg_group_coords <- direction
    length_vector_group <- sqrt(drop(avg_group_coords %*% avg_group_coords))

    # The line sometimes point into the "wrong" direction.
    # We can determine the cosine, if its negative we flip the x coords.

    cosangle <- 1

    if (!is.null(group)) {
        subgroup <- cent[group, ]

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
#' Calculates matrix of apl coordinates for random directions
#'
#' @inheritParams get_apl_cutoff
#' @param dims Number of dimensions to use for the analysis.
#'
#' @returns
#' List with permuted apl coordinates ("apl_perm") and, a list of saved ca
#' components ("saved_ca") that allow for quick recomputation of the CA results.
#'  For random_direction_cutoff this saved_ca is empty.
random_direction_cutoff <- function(
    caobj,
    dims = caobj@dims,
    apl_cutoff_reps = 300) {
    row_num <- nrow(caobj@prin_coords_cols)

    if (caobj@dims == 1 && !is.empty(caobj@dims)) {
        row_num <- 1
    }

    if (dims < caobj@dims) {
        caobj <- APL::subset_dims(caobj = caobj, dims = dims)
    }

    cols <- caobj@prin_coords_cols
    apl_perm <- data.frame(
        "x" = rep(0, row_num * apl_cutoff_reps),
        "y" = rep(0, row_num * apl_cutoff_reps)
    )


    for (k in seq(apl_cutoff_reps)) {
        # this picks random directions within the space of the original data.
        avg_group_coords <- stats::runif(
            n = dims,
            min = 0,
            max = stats::quantile(cols, 0.99)
        )

        length_vector_group <- sqrt(drop(avg_group_coords %*% avg_group_coords))
        length_vector_cols <- sqrt(rowSums(cols^2))

        colx <- drop(cols %*% avg_group_coords) / length_vector_group
        # pythagoras, y(r)=b²=c²-a²
        coly <- sqrt(length_vector_cols^2 - colx^2)

        colx[is.na(colx)] <- 0
        coly[is.na(coly)] <- 0

        idx <- ((1:row_num) + ((k - 1) * row_num))
        apl_perm[idx, ] <- cbind("x" = colx, "y" = coly)
    }
    return(apl_perm)
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
#' @inherit random_direction_cutoff return
#'
permutation_cutoff <- function(caobj,
                               counts,
                               group = caobj@group,
                               dims = caobj@dims,
                               apl_cutoff_reps = 10,
                               python = TRUE) {
    row_num <- nrow(counts)

    apl_perm <- data.frame(
        "x" = rep(0, row_num * apl_cutoff_reps),
        "y" = rep(0, row_num * apl_cutoff_reps)
    )


    if (caobj@dims == 1 && !is.empty(caobj@dims)) {
        row_num <- 1
    }

    margin <- 1
    pc <- 1
    cc <- TRUE
    cr <- FALSE

    for (k in seq(apl_cutoff_reps)) {
        # permute rows and rerun cacomp

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

        idx <- ((seq_len(row_num) + ((k - 1) * row_num)))

        apl_perm[idx, ] <- caobjp@apl_cols
    }

    return(apl_perm)
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
#'
#' @returns
#' The cutoff angle alpha in radians.
get_apl_cutoff <- function(caobj,
                           method = "random",
                           group = caobj@group,
                           counts = NULL,
                           quant = 0.99,
                           apl_cutoff_reps = 100) {
    if (method == "random") {
        if (is.null(apl_cutoff_reps)) {
            apl_cutoff_reps <- 100
        } else if (apl_cutoff_reps < 100) {
            warning("Number of repetitions should be set >=100.")
            apl_cutoff_reps <- 100
        }

        apl_perm <- random_direction_cutoff(
            caobj = caobj,
            dims = caobj@dims,
            apl_cutoff_reps = apl_cutoff_reps
        )
    } else if (method == "permutation") {
        if (is.null(apl_cutoff_reps)) {
            apl_cutoff_reps <- 5
        } else if (apl_cutoff_reps > 10) {
            message("Large number of repetitions might take a long time.")
        }

        apl_perm <- permutation_cutoff(
            caobj = caobj,
            counts = counts,
            group = group,
            dims = caobj@dims,
            apl_cutoff_reps = apl_cutoff_reps,
            python = TRUE
        )
    }

    # cotan between row and x axis
    apl_perm[, 3] <- apl_perm[, 1] / apl_perm[, 2]
    apl_perm[, 3][is.na(apl_perm[, 3])] <- 0

    cutoff_cotan <- stats::quantile(apl_perm[, 3], quant)

    # angle alpha is in radian.
    alpha <- atan(1 / cutoff_cotan)

    return(alpha)
}

# FIXME: DELETE ONLY USED BY get_apl_mergers
#' Calculates APL coordinates for points (columns)
#'  and directions of a cadir object.
#'  point into a precomputed APL space.
#' @param caobj A cacomp object.
#' @param cadir A cadir object.
#' @param apl_dir The direction in the original space for which to
#'  compute the APL coordinates.
#' @param group A vector of group indices.
#'  (Usually the cluster belonging to `apl_dir`)
#' @returns
#' A list with the columns (points) and directions projected into the APL space.
# apl_dir_coords <- function(cadir, caobj, apl_dir, group) {
#     model <- apl_model(
#         caobj = caobj,
#         direction = apl_dir,
#         group = group
#     )
#
#     apl_cols <- model(caobj@prin_coords_cols)
#     apl_dirs <- model(cadir@directions)
#     for (r in seq_len(nrow(apl_dirs))) {
#         sel <- match(
#             names(cadir@cell_clusters)[cadir@cell_clusters == rownames(apl_dirs)[r]],
#             rownames(caobj@prin_coords_cols)
#         )
#
#         if (length(sel) > 1) {
#             grp_mean <- colMeans(apl_cols[sel, ])
#         } else {
#             grp_mean <- apl_cols[sel, ]
#         }
#
#         if (sign(grp_mean[1]) != sign(apl_dirs[r, 1])) {
#             apl_dirs[r, ] <- c(-1, 1) * apl_dirs[r, ]
#         }
#     }
#
#     return(list("apl_cols" = apl_cols, "apl_dirs" = apl_dirs))
# }

# FIXME: DELETE. NOT NEEDED ANYMORE.
#' #' Checks if two clusters should be merged based on their angle in APL space.
#' #' @param cadir A cadir object.
#' #' @param caobj A cacomp object.
#' #' @param apl_quant The quantile to use for the cutoff.
#' #' @inheritParams get_apl_cutoff
#' #' @returns
#' #' A matrix with TRUE/FALSE values indicating if two clusters should be merged.
#' #' Importantly the function returns as soon as a candidate is found,
#' #'  the matrix stays the same size however.
#' get_apl_mergers <- function(cadir,
#'                             caobj,
#'                             apl_cutoff_reps = 100,
#'                             method = "random",
#'                             counts = NULL,
#'                             apl_quant = 0.99) {
#'     fun_args <- match.call()
#'     candidates <- matrix(FALSE,
#'         nrow = nrow(cadir@directions),
#'         ncol = nrow(cadir@directions)
#'     )
#'     dir_nms <- rownames(cadir@directions)
#'     for (d in seq_len(nrow(cadir@directions))) {
#'         grp_idx <- which(cadir@cell_clusters == dir_nms[d])
#'
#'         if (length(grp_idx) == 0) next
#'
#'         aplcds <- apl_dir_coords(
#'             cadir = cadir,
#'             caobj = caobj,
#'             apl_dir = cadir@directions[d, ],
#'             group = grp_idx
#'         )
#'
#'         cutoff_exists <- is_stored(cadir = cadir, fun_args = fun_args)
#'
#'         if (isTRUE(cutoff_exists)) {
#'             cutoff <- cadir@parameters$sa_cutoff
#'         } else {
#'             cutoff <- get_apl_cutoff(
#'                 caobj = caobj,
#'                 counts = counts,
#'                 method = method,
#'                 group = grp_idx,
#'                 apl_cutoff_reps = apl_cutoff_reps,
#'                 quant = apl_quant
#'             )
#'
#'             cadir@parameters$sa_cutoff <- cutoff
#'         }
#'
#'         cutoff <- rad_to_ang_sim(cutoff)
#'
#'         sim <- get_ang_sim(aplcds$apl_dirs[d, ], aplcds$apl_dirs)
#'         sim[1, d] <- 0
#'
#'         candidates[d, ] <- sim > cutoff
#'
#'         if (any(candidates[d, ])) {
#'             return(candidates)
#'         }
#'     }
#'
#'     return(candidates)
#' }

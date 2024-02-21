#TODO: saving cutoff results.

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
    reps = 300
) {
    row_num <- nrow(caobj@prin_coords_cols)

    if (caobj@dims == 1 && !is.empty(caobj@dims)) {
        row_num <- 1
    }

    if (dims < caobj@dims) {
        caobj <- APL::subset_dims(caobj = caobj, dims = dims)
    }

    cols <- caobj@prin_coords_cols
    apl_perm <- data.frame(
        "x" = rep(0, row_num * reps),
        "y" = rep(0, row_num * reps)
    )


    for (k in seq(reps)) {

        avg_group_coords <- runif(n = dims, min = 0, max = quantile(cols, 0.99))
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


# TODO: Add way to store permutation results.

#' Calculates permuted association plot coordinates
#'
#' @description
#' Calculates matrix of apl coordinates when permuting the original data.
#'
#' @inheritParams get_apl_cutoff
#' @param python If TRUE, use python for CA SVD.
#'
#' @inherit random_direction_cutoff return
#'
permutation_cutoff <- function(caobj,
                               counts,
                               group = caobj@group,
                               dims = caobj@dims,
                               reps = 10,
                               python = TRUE) {
    row_num <- nrow(counts)

    apl_perm <- data.frame(
        "x" = rep(0, row_num * reps),
        "y" = rep(0, row_num * reps)
    )


    if (caobj@dims == 1 && !is.empty(caobj@dims)) {
        row_num <- 1
    }

    margin <- 1
    pc <- 1
    cc <- TRUE
    cr <- FALSE

    for (k in seq(reps)) {
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
#' @param apl_cols A matrix of association plot coordinates.
#' @param method Method to use for computing the cutoff.
#' Either "random" or "permutation".
#' @param group A vector of group indices.
#' @param counts The original count matrix that was used to compute the caobj.
#' @param quant The quantile to use for the cutoff.
#' @param reps The number of repetitions to use.
#' Should be between 3-10 for permutation, and >=100 for random.
#'
#' @returns
#' The cutoff angle alpha in radians.
#'
get_apl_cutoff <- function(caobj,
                           method = "random",
                           group = caobj@group,
                           counts = NULL,
                           quant = 0.99,
                           reps = NULL) {

    if (method == "random") {

        if (is.null(reps)) {
            reps <- 100
        } else if (reps < 100) {
            warning("Number of repetitions should be set >=100.")
            reps <- 100
        }

        apl_perm <- random_direction_cutoff(
            caobj = caobj,
            dims = caobj@dims,
            reps = reps
        )

    } else if (method == "permutation") {
        if (is.null(reps)) {
            reps <- 5
        } else if (reps > 10) {
            message("Large number of repetitions might take a long time.")
        }

        apl_perm <- permutation_cutoff(
            caobj = caobj,
            mat = counts,
            group = group,
            dims = caobj@dims,
            reps = reps,
            store_perm = FALSE,
            python = TRUE
        )
    }

    # cotan between row and x axis
    apl_perm[, 3] <- apl_perm[, 1] / apl_perm[, 2]
    apl_perm[, 3][is.na(apl_perm[, 3])] <- 0

    cutoff_cotan <- quantile(apl_perm[, 3], quant)

    # angle alpha is in radian.
    alpha <- atan(1 / cutoff_cotan)

    return(alpha)
}


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
apl_dir_coords <- function(cadir, caobj, apl_dir, group) {

    model <- apl_model(
        caobj = caobj,
        direction = apl_dir,
        group = group
    )

    apl_cols <- model(caobj@prin_coords_cols)
    apl_dirs <- model(cadir@directions)

    for (r in seq_len(nrow(apl_dirs))) {

        sel <- match(
            names(cadir@cell_clusters)[cadir@cell_clusters == r],
            rownames(caobj@prin_coords_cols)
        )

        if (length(sel) > 1) {
            grp_mean <- colMeans(apl_cols[sel, ])
        } else {
            grp_mean <- apl_cols[sel, ]
        }

        if (sign(grp_mean[1]) != sign(apl_dirs[r, 1])) {
            apl_dirs[r, ] <- c(-1, 1) * apl_dirs[r, ]
        }

    }

    return(list("apl_cols" = apl_cols, "apl_dirs" = apl_dirs))

}

#' Checks if two clusters should be merged based on their angle in APL space.
#' @param cadir A cadir object.
#' @param caobj A cacomp object.
#' @param apl_quant The quantile to use for the cutoff.
#' @inheritParams get_apl_cutoff
#' @inheritParams apl_dir_cutoff
#' @returns
#' A matrix with TRUE/FALSE values indicating if two clusters should be merged.
#' Importantly the function returns as soon as a candidate is found,
#'  the matrix stays the same size however.
get_apl_mergers <- function(cadir,
                            caobj,
                            reps = 100,
                            method = "random",
                            counts = NULL,
                            apl_quant = 0.99) {

    candidates <- matrix(FALSE,
                         nrow = nrow(cadir@directions),
                         ncol = nrow(cadir@directions))

    for (d in seq_len(nrow(cadir@directions))) {

        grp_idx <- which(f2n(cadir@cell_clusters) == d)

        if (length(grp_idx) == 0) next

        aplcds <- apl_dir_coords(cadir = cadir,
                                 caobj = caobj,
                                 apl_dir = cadir@directions[d, ],
                                 group = grp_idx)


        cutoff <- get_apl_cutoff(caobj = caobj,
                                 counts = counts,
                                 method = method,
                                 group = grp_idx,
                                 reps = reps,
                                 quant = apl_quant)

        cutoff <- rad_to_ang_sim(cutoff)

        sim <- get_ang_sim(aplcds$apl_dirs[d, ], aplcds$apl_dirs)
        sim[1, d] <- 0

        candidates[d, ] <- sim > cutoff

        if (any(candidates[d, ])) {
            return(candidates)
        }
    }

    return(candidates)
}



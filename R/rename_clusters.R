#' Rename clusters and directions so that they match.
#' @param cadir A cadir object.
#' @details
#' rename_clusters takes a cadir object and renames the clusters and directions
#' such that they are again coherent (e.g. from 1:5).
#' @returns
#' A cadir object with renamed clusters and directions.
rename_clusters  <- function(cadir) {

    cc_nms <- levels(cadir@cell_clusters)
    dir_nms <- rownames(cadir@directions)

    dict <- cadir@dict

    # Remove directions without any CELLS clustered.
    # There might be directions that only have genes,
    # but we remove it nevertheless.
    have_smpls <- which(dir_nms %in% cc_nms)
    cadir@directions <- cadir@directions[have_smpls, , drop = FALSE]
    dir_nms <- rownames(cadir@directions)

    # NOTE: We first fix the names of the directions.
    # Then we can use those to fix the cell/gene cluster names.
    # Also, we want to fix the dict around this time too.

    # 1) Rename the directions according to the overall ordering
    # But leave none-standard cluster names alone.

    # Subset dict to entries for which directions exist.
    # Get the positions of dir cluster names in dict.
    dict <- dict[names(dict) %in% dir_nms]
    dict <- dict[!duplicated(names(dict))]
    dict_pos <- match(dir_nms, names(dict))

    # Find directions that are standard naming.
    is_std_dir <- is_std_name(dir_nms)

    # Change to name based on position in matrix
    old_dir_nms <- dir_nms
    dir_nms[is_std_dir] <- cl2nm(which(is_std_dir))
    rownames(cadir@directions) <- dir_nms

    # 2) Change cluster names corresponding to dir in dict.
    names(dict)[dict_pos] <- dir_nms
    dict[dir_nms] <- seq_len(length(dir_nms))
    cadir@dict <- dict

    # 3) Change cell cluster names
    new_cls_nms <- dir_nms[match(cadir@cell_clusters, old_dir_nms)]
    cadir@cell_clusters <- stats::setNames(
        factor(new_cls_nms, levels = dir_nms),
        names(cadir@cell_clusters)
    )

    # 4) Change gene cluster names
    if (!is.empty(cadir@gene_clusters)) {
        new_cls_nms <- dir_nms[match(cadir@gene_clusters, old_dir_nms)]
        cadir@gene_clusters <- stats::setNames(
            factor(new_cls_nms, levels = dir_nms),
            names(cadir@gene_clusters)
        )
    }

    # 5) Change distances
    if (!is.empty(cadir@distances)) {
        # cadir@distances <- cadir@distances[, have_smpls, drop = FALSE]
        # colnames(cadir@distances) <- dir_nms
    }

    # Double Check that the dict is correct!
    if (isFALSE(.check_dict(cadir))) {
        cadir <- .correct_dict(cadir)
    }

    # FIXME: Double check that combined cell & gene levels are not a problem!
    cadir@cell_clusters <- droplevels(cadir@cell_clusters)
    cadir@gene_clusters <- droplevels(cadir@gene_clusters)
    comb_lvl <- unique(c(levels(cadir@cell_clusters), levels(cadir@gene_clusters)))
    levels(cadir@cell_clusters) <- comb_lvl
    levels(cadir@gene_clusters) <- comb_lvl

    stopifnot(!any(is.na(cadir@cell_clusters)))
    stopifnot(!any(is.na(cadir@gene_clusters)))

    return(cadir)
}

#' Runs an array of quick tests to check if the @dict slot is coherant.
#' @param cadir A CAdir object with a @dict slot.
#' @returns
#' TRUE if everything is correct in the dict.
.check_dict <- function(cadir) {

    ord <- match(rownames(cadir@directions), names(cadir@dict))
    correct_order <- all(
        unlist(cadir@dict[ord]) == seq_len(nrow(cadir@directions))
    )

    all_ccs <- all(unique(cadir@cell_clusters %in% names(cadir@dict)))
    all_gcs <- all(unique(cadir@gene_clusters %in% names(cadir@dict)))

    no_extra <- all(names(cadir@dict) %in% rownames(cadir@directions))
    no_missing <- all(rownames(cadir@directions) %in% names(cadir@dict))
    no_dupl <- anyDuplicated(names(cadir@dict))

    all_correct <- correct_order &&
        all_ccs &&
        all_gcs &&
        no_extra &&
        no_missing &&
        no_dupl

    return(all_correct)
}

#' Attempts to correct a malformed dict slot.
#' @param cadir A CAdir object with a @dict slot.
#' @returns
#' A CAdir object with corrected dict slot.
.correct_dict <- function(cadir) {

    dict <- cadir@dict

    # Make unique
    dict <- dict[!duplicated(names(dict))]

    # Check for dict entries that dont exist
    to_keep <- which(names(dict) %in% rownames(cadir@directions))
    dict <- dict[to_keep]

    # Check for dict entries that are missing
    to_add <- which(!rownames(cadir@directions) %in% names(dict))
    names(to_add) <- rownames(cadir@directions)[to_add]

    dict <- c(
        dict,
        as.list(to_add)
    )

    # Order the dict.
    ord <- match(rownames(cadir@directions), names(cadir@dict))
    dict <- dict[ord]

    # Some extra checks
    stopifnot(all(unique(cadir@cell_clusters %in% names(dict))))
    stopifnot(all(unique(cadir@gene_clusters %in% names(dict))))

    cadir@dict <- dict
    return(cadir)

}

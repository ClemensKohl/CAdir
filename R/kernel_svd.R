# https://bobrupakroy.medium.com/what-is-kernel-pca-using-r-python-4864c2471e62

# FIXME: WIP
kernel_ca <- function(obj,
                      dims = 2,
                      top = nrow(obj),
                      residuals = "pearson",
                      clip = FALSE,
                      cutoff = NULL,
                      kernel = "rbfdot",
                      kpar = list(sigma = 0.1)) {

    obj <- APL:::rm_zeros(obj)

    if (!is.null(top) && top < nrow(obj)) {

        obj <- APL:::var_rows(mat = obj,
                              top = top,
                              residuals = residuals,
                              clip = clip,
                              cutoff = cutoff)
    }

    top <- nrow(obj)

    if (top > nrow(obj)) {
        warning("\nParameter top is >nrow(obj) and therefore ignored.")
    } else if (is.null(top) || top == nrow(obj)) {
        # do nothing. just here to allow for else statement.
    } else {
        warning("\nUnusual input for top, argument ignored.")
    }

    res <- APL:::calc_residuals(mat = obj,
                                residuals = residuals,
                                clip = clip,
                                cutoff = cutoff)
    s <- t(res$S)

    kpca <- kernlab::kpca(x = s,
                          kernel = kernel,
                          features = dims)

    # S = UDV^T
    # cxg = [cxd][dxd][dxg]

    v <- kpca@pcv
    e <- kpca@eig
    d <- sqrt(e) # D <- sqrt(E)*(nrow(S) - 1)
    u <- t(diag(1 / d) %*% t(v) %*% kpca@xmatrix) # TODO: Check if correct


    names(d) <- paste0("Dim", seq_len(length(d)))
    dimnames(v) <- list(rownames(s), paste0("Dim", seq_len(ncol(v))))
    dimnames(u) <- list(colnames(s), paste0("Dim", seq_len(ncol(u))))


    caobj <- do.call(
        APL::new_cacomp,
        list(
            "U" = u,
            "V" = v,
            "D" = d,
            "row_masses" = res$rowm,
            "col_masses" = res$colm,
            "top_rows" = top,
            "tot_inertia" = res$tot,
            "dims" = dims
        )
    )

    caobj <- APL::ca_coords(caobj)

    return(caobj)
}

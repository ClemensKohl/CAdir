# https://bobrupakroy.medium.com/what-is-kernel-pca-using-r-python-4864c2471e62

# setMethod(
#     "kpca", signature(x = "matrix"),
#     function(x, kernel = "rbfdot", kpar = list(sigma = 0.1), features = 0, th = 1e-4, na.action = na.omit, ...) {
#         x <- na.action(x)
#         x <- as.matrix(x)
#         m <- nrow(x)
#         ret <- new("kpca")
#         if (!is(kernel, "kernel")) {
#             if (is(kernel, "function")) kernel <- deparse(substitute(kernel))
#             kernel <- do.call(kernel, kpar)
#         }
#         if (!is(kernel, "kernel")) stop("kernel must inherit from class `kernel'")
#
#         km <- kernelMatrix(kernel, x)
#
#         ## center kernel matrix
#         kc <- t(t(km - colSums(km) / m) - rowSums(km) / m) + sum(km) / m^2
#
#         ## compute eigenvectors
#         res <- eigen(kc / m, symmetric = TRUE)
#
#         if (features == 0) {
#             features <- sum(res$values > th)
#         } else if (res$values[features] < th) {
#             warning(paste("eigenvalues of the kernel matrix are below threshold!"))
#         }
#
#         pcv(ret) <- t(t(res$vectors[, 1:features]) / sqrt(res$values[1:features]))
#         eig(ret) <- res$values[1:features]
#         names(eig(ret)) <- paste("Comp.", 1:features, sep = "")
#         rotated(ret) <- kc %*% pcv(ret)
#         kcall(ret) <- match.call()
#         kernelf(ret) <- kernel
#         xmatrix(ret) <- x
#         return(ret)
#     }
# )
loadNamespace("kernlab")
setClass("rbfkernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))
rbfdot <- function(sigma = 1) {
    rval <- function(x, y = NULL) {
        if (!is(x, "vector")) stop("x must be a vector")
        if (!is(y, "vector") && !is.null(y)) stop("y must a vector")
        if (is(x, "vector") && is.null(y)) {
            return(1)
        }
        if (is(x, "vector") && is(y, "vector")) {
            if (!length(x) == length(y)) {
                stop("number of dimension must be the same on both data points")
            }
            return(exp(sigma * (2 * crossprod(x, y) - crossprod(x) - crossprod(y))))
            # sigma/2 or sigma ??
        }
    }
    return(new("rbfkernel", .Data = rval, kpar = list(sigma = sigma)))
}

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
    # s <- t(res$S)

    # FIXME: The function centers the data, we dont want that.
    # kpca <- kernlab::kpca(x = s,
    #                       kernel = kernel,
    #                       features = dims,
    #                       kpar = kpar)

    s <- as.matrix(t(res$S))
    kern <- do.call(kernel, kpar)
    km <- kernlab::kernelMatrix(kernel = kern, x = s)

    ## center kernel matrix
    # kc <- t(t(km - colSums(km) / m) - rowSums(km) / m) + sum(km) / m^2

    ## compute eigenvectors
    m <- nrow(s)
    kpca <- eigen(km / m, symmetric = TRUE)

    v <- t(t(kpca$vectors[, 1:dims]) / sqrt(kpca$values[1:dims]))
    e <- kpca$values[1:dims]
    d <- sqrt(e) # D <- sqrt(E)*(nrow(S) - 1)
    names(e) <- paste("Comp.", 1:dims, sep = "")
    prin_coords_cols <- km %*% v

    # NOTE: We have to provide a matrix U, so we make a dummy vector.
    u <- matrix(NA, nrow = ncol(s), ncol = dims)

    # eigenvalues in a decreasing order
    # svd_res <- irlba::irlba(A = km, nv = dims, smallest = FALSE)
    # svd_res <- svd_res[1:3]
    # names(svd_res)[1:3] <- c("D", "U", "V")
    # svd_res$D <- as.vector(svd_res$D)

    # S = UDV^T
    # cxg = [cxd][dxd][dxg]

    # v <- kpca@pcv
    # e <- kpca@eig
    # d <- sqrt(e) # D <- sqrt(E)*(nrow(S) - 1)
    # u <- t(diag(1 / d) %*% t(v) %*% kpca@xmatrix) # NOTE: This cannot be done! PCA is in kernel space!

    names(d) <- paste0("Dim", seq_len(length(d)))
    dimnames(v) <- list(rownames(s), paste0("Dim", seq_len(ncol(v))))
    dimnames(u) <- list(colnames(s), paste0("Dim", seq_len(ncol(u))))


    caobj <- do.call(
        APL::new_cacomp,
        list(
            "U" = u,
            "V" = v,
            "D" = d,
            "prin_coords_cols" = prin_coords_cols,
            "row_masses" = res$rowm,
            "col_masses" = res$colm,
            "top_rows" = top,
            "tot_inertia" = res$tot,
            "dims" = dims
        )
    )

    # caobj <- APL::ca_coords(caobj) # NOTE: Again, we cannot compute this, as we are working in kernel space!

    return(caobj)
}

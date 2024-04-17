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

    data <- APL:::rm_zeros(obj)

    if (!is.null(top) && top < nrow(obj)) {

        obj <- APL:::var_rows(mat = obj,
                        top = top,
                        residuals = residuals,
                        clip = clip,
                        cutoff = cutoff)
        toptmp <- top
    }

    top <- toptmp <- nrow(obj)

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
    S <- t(res$S)

    kpca <- kernlab::kpca(x = S,
                          kernel = kernel,
                          features = dims)

    # X = UDV^T
    # cxg = [cxd][dxd][dxg]

    V <- kpca@pcv
    E <- kpca@eig
    D <- sqrt(E) # D <- sqrt(E)*(nrow(S) - 1)
    U <- t(diag(1/D) %*% t(V) %*% kpca@xmatrix) # TODO: Check if correct


    names(D) <- paste0("Dim", seq_len(length(D)))
    dimnames(V) <- list(rownames(S), paste0("Dim", seq_len(ncol(V))))
    dimnames(U) <- list(colnames(S), paste0("Dim", seq_len(ncol(U))))


    caobj <- do.call(APL::new_cacomp,
                     list(
                          "U" = U,
                          "V" = V,
                          "D" = D,
                          "row_masses" = res$rowm,
                          "col_masses" = res$colm,
                          "top_rows" = top,
                          "tot_inertia" = res$tot,
                          "dims" = dims
                          ))

    caobj <- APL::ca_coords(caobj)

    return(caobj)
}


## TESTING
library(SingleCellExperiment)
library(kernlab)
library(RUtils)

sce <- RUtils:::get_zeisel_brain()

kca <- kernel_ca(obj = logcounts(sce),
                 dims = 30,
                 top = 2000,
                 residuals = "pearson",
                 kernel = "rbfdot",
                 kpar = list(sigma = 0.3))

kca <- kernel_ca(obj = logcounts(sce),
                 dims = 30,
                 top = 2000,
                 residuals = "pearson",
                 kernel = "vanilladot")

ca <- APL::cacomp(obj = logcounts(sce),
                   dims = 30,
                   top = 2000,
                   residuals = "pearson")

cabic <- dirclust_splitmerge(caobj = kca,
                             k = 7,
                             cutoff = 85,
                             qcutoff = 0.8,
                             method = "random",
                             counts = NULL,
                             apl_quant = 0.99,
                             min_cells = 5,
                             epochs = 15,
                             reps = 5,
                             make_plots = FALSE)

APL::ca_biplot(kca)


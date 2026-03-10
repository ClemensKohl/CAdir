# Reference implementations of the old code used to verify refactored versions.

ref_random_cotan <- function(caobj, dims, apl_cutoff_reps, axis = "cols") {
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
  apl_perm <- data.frame(
    "x" = rep(0, row_num * apl_cutoff_reps),
    "y" = rep(0, row_num * apl_cutoff_reps)
  )

  for (k in seq(apl_cutoff_reps)) {
    avg_group_coords <- stats::runif(
      n = dims,
      min = 0,
      max = stats::quantile(cols, 0.99)
    )

    length_vector_group <- sqrt(drop(avg_group_coords %*% avg_group_coords))
    length_vector_coords <- sqrt(rowSums(ca_coords^2))

    coord_x <- drop(ca_coords %*% avg_group_coords) / length_vector_group
    coord_y <- sqrt(length_vector_coords^2 - coord_x^2)

    coord_x[is.na(coord_x)] <- 0
    coord_y[is.na(coord_y)] <- 0

    idx <- ((1:row_num) + ((k - 1) * row_num))
    apl_perm[idx, ] <- cbind("x" = coord_x, "y" = coord_y)
  }

  cotan <- apl_perm[, 1] / apl_perm[, 2]
  cotan[is.na(cotan)] <- 0
  return(cotan)
}


ref_permutation_cotan <- function(
  caobj,
  counts,
  group,
  dims,
  apl_cutoff_reps,
  axis = "cols"
) {
  if (axis == "cols") {
    row_num <- nrow(caobj@V)
  } else if (axis == "rows") {
    row_num <- nrow(caobj@U)
  }

  apl_perm <- data.frame(
    "x" = rep(0, row_num * apl_cutoff_reps),
    "y" = rep(0, row_num * apl_cutoff_reps)
  )

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
        python = FALSE,
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
    if (axis == "cols") {
      apl_perm[idx, ] <- caobjp@apl_cols
    } else if (axis == "rows") {
      apl_perm[idx, ] <- caobjp@apl_rows
    }
  }

  cotan <- apl_perm[, 1] / apl_perm[, 2]
  cotan[is.na(cotan)] <- 0
  return(cotan)
}


make_test_caobj <- function(seed = 5, n_genes = 50, n_cells = 20, dims = 5) {
  set.seed(seed)
  counts <- matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes,
    ncol = n_cells
  )
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- paste0("cell", seq_len(n_cells))
  caobj <- APL::cacomp(
    obj = counts,
    dims = dims,
    top = n_genes,
    residuals = "pearson",
    princ_coords = 3
  )
  list(caobj = caobj, counts = counts)
}


test_that("random_direction_cutoff returns same cotangents as reference", {
  ca <- make_test_caobj()
  caobj <- ca$caobj
  dims <- caobj@dims
  reps <- 20

  set.seed(42)
  ref <- ref_random_cotan(caobj, dims = dims, apl_cutoff_reps = reps)

  set.seed(42)
  new <- random_direction_cutoff(caobj, dims = dims, apl_cutoff_reps = reps)

  expect_equal(new, ref)
})


test_that("permutation_cutoff returns same cotangents as reference", {
  ca <- make_test_caobj()
  caobj <- ca$caobj
  counts <- ca$counts
  group <- seq_len(10)
  dims <- caobj@dims
  reps <- 2

  set.seed(42)
  ref <- ref_permutation_cotan(
    caobj = caobj,
    counts = counts,
    group = group,
    dims = dims,
    apl_cutoff_reps = reps
  )

  set.seed(42)
  new <- permutation_cutoff(
    caobj = caobj,
    counts = counts,
    group = group,
    dims = dims,
    apl_cutoff_reps = reps,
    python = FALSE
  )

  expect_equal(new, ref)
})


test_that("get_apl_cutoff random method returns finite positive alpha", {
  ca <- make_test_caobj()

  set.seed(1)
  alpha <- get_apl_cutoff(ca$caobj, method = "random", apl_cutoff_reps = 100)

  expect_true(is.finite(alpha))
  expect_gt(alpha, 0)
})

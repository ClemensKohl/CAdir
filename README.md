# Clustering by directions in CA space

This package implements a clustering algorithms that determines clusters iteratively by their directions in CA space.
Unlike most other algorithms it does not require prior knowledge of the number of clusters in the data, but can instead infer them during clustering.

# TODO
- [ ] Add documentation to missing functions.
- [ ] Reduce dependencies and insure that all packages are available on CRAN/Bioconductor. -> ggsankey!
- [x] "Prettify" plotting functions. -> Colors for APL plots!
- [ ] Write vignette.
- [x] Add Salpha scoring to rank the genes.
- [ ] Add kernel CA as a non-linear alternative. (Prototype exists on branch kernel_svd)
- [ ] Enable labeling of genes in cluster_apl

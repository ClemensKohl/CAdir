# Clustering by directions in CA space

This package implements a clustering algorithms that determines clusters iteratively by their directions in CA space.
Unlike most other algorithms it does not require prior knowledge of the number of clusters in the data, but can instead infer them during clustering.

The package is fully functioning, but it is still work in progress and edge cases etc. still have to be ironed out.

# TODO

- [x] "Prettify" plotting functions. -> Colors for APL plots!
- [x] Add Salpha scoring to rank the genes.
- [x] Enable labeling of genes in cluster_apl
- [x] **Split anti-correlated clusters!**
- [x] **Fix naming of directions. Should be coherent with cluster names at all times!**
- [ ] Make a verbose toggle.
- [ ] Add missing documentation to functions.
- [ ] Add kernel CA as a non-linear alternative. (Prototype exists on branch kernel_svd)
- [ ] Write vignette.
- [ ] Reduce dependencies and insure that all packages are available on CRAN/Bioconductor. -> ggsankey!


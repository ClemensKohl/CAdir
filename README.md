# Clustering by directions in CA space

This package implements a clustering algorithms that determines clusters iteratively by their directions in CA space.
Unlike most other algorithms it does not require prior knowledge of the number of clusters in the data, but can instead infer them during clustering.
The package can be installed through GitHub:

``` r
devtools::install_github("ClemensKohl/CAdir")
```

Note that currently you have to also install CAbiNet from GitHub for the package to function.

# Quick start

Download example data and perform basic preprocessing:

``` r
library(scRNAseq)
library(scran)
library(scater)
library(scuttle)

set.seed(2358)

sce <- scRNAseq::ZeiselBrainData()
clust <- scran::quickCluster(sce)
sce <- scran::computeSumFactors(sce, cluster = clust, min.mean = 0.1)
sce <- scuttle::logNormCounts(sce)
dec <- scran::modelGeneVar(sce)
top_genes <- scran::getTopHVGs(dec, prop = 0.8)
sce <- sce[top_genes, ]
```

## CAdir

``` r
cak <- dirclust_splitmerge(
    caobj = ca,
    k = 2,
    cutoff = NULL,
    method = "random",
    apl_quant = 0.99,
    counts = NULL,
    min_cells = 20,
    reps = 4,
    make_plots = TRUE,
    apl_cutoff_reps = 100
)

```

``` r
cadir <- rank_genes(cadir = cak, caobj = ca)
top <- top_genes(cadir)

anncak <- annotate_biclustering(obj = cak,
                                universe = rownames(sce_sub),
                                org = "mm")
```

``` r
cluster_apl(cadir = cak,
            caobj = ca,
            direction = cak@directions[4,],
            group = which(cak@cell_clusters == "cluster_4"),
            cluster = "cluster_4",
            show_genes = TRUE,
            show_cells = TRUE,
            show_lines = FALSE,
            highlight_cluster = FALSE)
```

``` r
sm_plot(cadir = cak,
        caobj = ca,
        rm_redund = TRUE,
        keep_end = TRUE,
        highlight_cluster = F,
        show_genes = FALSE,
        show_cells = TRUE,
        annotate_clusters = T,
        org = "hs")
```

Verbosity of the messages can be controlled with rlang.
To turn all messages off:

``` r
rlang::local_options(mypackage.verbose = "quiet")
```

 To turn them back on:

``` r
rlang::local_options(mypackage.verbose = "verbose")
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_plots.R
\name{sm_plot}
\alias{sm_plot}
\title{Plots the graph of the clustering splits and merges
and overlays APL plots over the nodes.}
\usage{
sm_plot(
  cadir,
  caobj,
  rm_redund = TRUE,
  show_cells = TRUE,
  show_genes = FALSE,
  highlight_cluster = FALSE,
  annotate_clusters = FALSE,
  org = "mm",
  keep_end = TRUE,
  inlet_side = 0.1,
  title_size = 10,
  show_axis = FALSE,
  n_wrap = Inf,
  layout = "tree",
  rotate = TRUE
)
}
\arguments{
\item{cadir}{A cadir object for which to compute the APL}

\item{caobj}{A cacomp object.}

\item{rm_redund}{If TRUE, only shows an iteration if something changes.}

\item{show_cells}{If TRUE, points (cells) are plotted.}

\item{show_genes}{If TRUE, genes are plotted.}

\item{highlight_cluster}{If TRUE, highlights the points in belonging to
\code{cluster}, unless colour_by_group is TRUE.}

\item{annotate_clusters}{If TRUE uses automatic cell type annotation via
\code{annotate_biclustering} to annotate the clusters in each level.}

\item{org}{Organisme to use for the automatic annotation.
Either "mm" for mouse or "hs" for human.}

\item{keep_end}{Adds the last iteration, even if no change to previous.}

\item{inlet_side}{Length of the inlet side. The inlet is a square centered
around the cluster node.}

\item{title_size}{Font size of the title.}

\item{show_axis}{Whether or not to show axis markings.}

\item{n_wrap}{Character length around which to wrap inlet titles.}

\item{layout}{Graph layout.
Currently only supports \code{"dendrogram"} or \code{"tree"}.}

\item{rotate}{TRUE/FALSE. Whether to rotate the tree by 90 degrees.}
}
\value{
A ggplot object showing the split-merge graph and APL plots for each cluster.
}
\description{
Plots the graph of the clustering splits and merges
and overlays APL plots over the nodes.
}

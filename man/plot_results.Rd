% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_results}
\alias{plot_results}
\title{Plot for the clustering results wich shows the
relationship between the clusters. In the diagonal
an APL plot for the respective cluster is shown.}
\usage{
plot_results(
  cadir,
  caobj,
  highlight_cluster = TRUE,
  show_cells = TRUE,
  show_genes = FALSE,
  title_prefix = "Cluster: ",
  gsub_title = "",
  return_list = FALSE,
  ...
)
}
\arguments{
\item{cadir}{A cadir object with valid cell clustering results.}

\item{caobj}{A cacomp object.}

\item{highlight_cluster}{If TRUE, highlights the points in belonging to
\code{cluster}, unless colour_by_group is TRUE.}

\item{show_cells}{If TRUE, points (cells) are plotted.}

\item{show_genes}{If TRUE, genes are plotted.}

\item{title_prefix}{Prefix to print before cluster name.}

\item{gsub_title}{Character string to substitute with " " in plot titles.}

\item{return_list}{If TRUE, instead of a final, arranged panel the plots
are returned as a list.}

\item{...}{Further arguments forwarded to cluster_apl().}
}
\value{
A plot that summarizes the cell clustering results and
how the clusters relate to each other.
}
\description{
Plot for the clustering results wich shows the
relationship between the clusters. In the diagonal
an APL plot for the respective cluster is shown.
}

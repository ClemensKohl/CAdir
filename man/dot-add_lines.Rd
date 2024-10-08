% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_apl.R
\name{.add_lines}
\alias{.add_lines}
\title{Adds lines of other clusters to the APL.}
\usage{
.add_lines(ggplt, apl_dir, highlight_cluster = FALSE)
}
\arguments{
\item{ggplt}{More or less finished ggplot object.}

\item{apl_dir}{The cluster directions in APL coordinates.}

\item{highlight_cluster}{If TRUE, highlights the points in belonging to
\code{cluster}, unless colour_by_group is TRUE.}
}
\value{
ggplot with lines added.
}
\description{
Adds lines of other clusters to the APL.
}

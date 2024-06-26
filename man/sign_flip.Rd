% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dirclust.R
\name{sign_flip}
\alias{sign_flip}
\title{Determine sign for SVD singular vectors.}
\usage{
sign_flip(points, line)
}
\arguments{
\item{points}{Points/Cells for which to check the directionality of the line.}

\item{line}{Line/direction whose directionality needs to be determined.}
}
\value{
TRUE if the sign of the direction needs to be flipped, FALSE otherwise.
}
\description{
Determine sign for SVD singular vectors.
}
\references{
Bro, Rasmus, Acar, Evrim, and Kolda, Tamara Gibson.
Resolving the sign ambiguity in the singular value decomposition.
United States: N. p., 2007. Web. doi:10.2172/920802.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{Sgene2Pguides_hash}
\alias{Sgene2Pguides_hash}
\title{Gene hashing}
\usage{
Sgene2Pguides_hash(guide2gene, cores = 1)
}
\arguments{
\item{guide2gene}{derived from Input object}

\item{cores}{number of cores to use (default 1)}
}
\value{
Gene hash/list
}
\description{
Gene hashing
}
\examples{
\dontrun{
#' data("Input", package = "gemini")
Sgene2Pguides_hash(Input$guide.pair.annot[1:10,])
}
}

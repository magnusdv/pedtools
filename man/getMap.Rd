% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMap.R
\name{getMap}
\alias{getMap}
\alias{setMap}
\title{Tabulate marker positions}
\usage{
getMap(x, markers = NULL, na.action = 0, verbose = TRUE)

setMap(x, map, matchNames = NA, ...)
}
\arguments{
\item{x}{An object of class \code{ped} or a list of such.}

\item{markers}{A vector of names or indices referring to markers attached to
\code{x}. By default, all markers are included.}

\item{na.action}{Either 0 (default), 1 or 2. (See Details.)}

\item{verbose}{A logical.}

\item{map}{Either a data frame or the path to a map file.}

\item{matchNames}{A logical; if TRUE, pre-existing marker names of \code{x} will
be used to assign chromosome labels and positions from \code{map}.}

\item{...}{Further arguments passed to \code{read.table()}.}
}
\value{
\code{getMap()} returns a data frame with columns \code{CHROM}, \code{MARKER} and
\code{MB}.

\code{setMap()} returns \code{x} with modified marker attributes.
}
\description{
Return a map of the markers attached to a pedigree.
}
\details{
The \code{na.action} argument controls how missing values are dealt with:
\itemize{
\item \code{na.action} = 0: Return map unmodified
\item \code{na.action} = 1: Replace missing values with dummy values.
\item \code{na.action} = 2: Remove markers with missing data.
}

In \code{setMap()}, the \code{map} argument should be a data frame (or file) with the
following columns in order: 1) chromosome, 2) marker name, 3) position in megabases.
Column names are ignored, as are any columns after the first three.
}
\examples{
x = singleton(1)
m1 = marker(x, chrom = 1, posMb = 10, name = "m1")
m2 = marker(x, chrom = 1, posMb = 11)
m3 = marker(x, chrom = 1)
x = setMarkers(x, list(m1, m2, m3))

# Compare effect of `na.action`
getMap(x, na.action = 0)
getMap(x, na.action = 1)
getMap(x, na.action = 2)

# Getting and setting map are inverses
y = setMap(x, getMap(x))
identical(x,y)

}

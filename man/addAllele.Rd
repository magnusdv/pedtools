% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marker_utils.R
\name{addAllele}
\alias{addAllele}
\title{Add allele}
\usage{
addAllele(x, marker, allele, freq = 0.001, adjust = c("previous", "all"))
}
\arguments{
\item{x}{A \code{ped} object or a list of such, or a frequency database (list of
numeric vectors).}

\item{marker}{The name or index of a marker attached to \code{x}.}

\item{allele}{The name of the new allele.}

\item{freq}{The frequency of the new allele, by default 0.001.}

\item{adjust}{Either "previous" or "all", indicating how the frequencies
should be adjusted so that they sum to 1. If "previous" (default), the
frequencies of the original alleles are multiplied with \code{1 - freq}. If
"all", scaling is performed after adding the new allele, i.e., dividing all
frequencies by \code{1 + freq}.}
}
\value{
A copy of \code{x} with modified marker attributes.
}
\description{
Extends the allele set of a marker attached to a pedigree, by adding a single
allele.
}
\examples{

## Ped input
x = nuclearPed() |>
  addMarker(geno = c(NA, NA, "b/c"), afreq = c(b = 0.5, c = 0.5))

y = addAllele(x, marker = 1, allele = "a")
afreq(y, 1)

z = addAllele(y, marker = 1, allele = "d", freq = 0.1, adjust = "all")
afreq(z, 1)


## Database input
db = list(M1 = c(a = .2, b = .3, c = .5),
          M2 = c("7" = .9, "8.3" = .1))
addAllele(db, marker = "M2", allele = "8")

}

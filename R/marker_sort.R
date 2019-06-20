#' Sort the alleles in each genotype
#'
#' Ensure that all genotypes are sorted internally. For example, if a marker
#' attached to `x` has alleles 1 and 2, then running this function will replace
#' all genotypes "2/1" by "1/2".
#'
#' @param x A `ped` object or a list of such
#'
#' @return An object identical to `x` except that the all genotypes are sorted.
#'
#' @examples
#' x = singleton(1)
#'
#' # Various markers with misordered genotypes
#' m1 = marker(x, `1` = 2:1)
#' m2 = marker(x, `1` = c('b','a'))
#' m3 = marker(x, `1` = c("100.3", "99.1"))
#' x = setMarkers(x, list(m1, m2, m3))
#' x
#'
#' # Sort all genotypes
#' y = sortGenotypes(x)
#' y
#'
#' # Also works when input is a list of peds
#' sortGenotypes(list(x, x))
#'
#' @export
sortGenotypes = function(x) {
  if(is.pedList(x))
    return(lapply(x, sortGenotypes))

  if(!hasMarkers(x)) return(x)

  # Need integer allele representation for comparisons.
  # Hence as.matrix.ped() instead of getAlleles().
  pedm = as.matrix(x)

  # Index of allele-1 columns (not hard-coding number of ped columns)
  a1idx = seq(ncol(pedm) - 2*nMarkers(x) + 1, ncol(pedm), by = 2)
  a2idx = a1idx + 1

  a1 = pedm[, a1idx]
  a2 = pedm[, a2idx]
  swap = a1 > a2

  pedm[, a1idx][swap] = a2[swap]
  pedm[, a2idx][swap] = a1[swap]

  restorePed(pedm)
}

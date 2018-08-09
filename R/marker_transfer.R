#' Transfer marker data
#'
#' Transfer marker data between pedigrees. Any markers attached to the
#' target are overwritten.
#'
#' @param from a `ped` or `singleton` object, or a
#'   list of such objects.
#' @param to a `ped` or `singleton` object, or a
#'   list of such objects.
#' @return A `ped` object (or a list of such) similar to `to`, but
#'   where all individuals also present in `from` have marker genotypes
#'   copied over.  Any previous marker data is erased.
#'
#' @examples
#'
#' x = nuclearPed(fa = "father", mo = "mother", children = "boy")
#' m = marker(x, father = 1:2, boy = 1:2)
#' x = setMarkers(x, m)
#'
#' y = list(singleton("father"), nuclearPed(children = "boy"))
#' y = transferMarkers(x, y)
#' y
#' stopifnot(genotype(y[[1]], 1, "father") == 1:2, genotype(y[[2]], 1, "boy") == 1:2)
#'
#' @export
transferMarkers = function(from, to) {
  if (is.ped(from) && is.ped(to))
    return(.transferMarkersSimple(from, to))
  if (is.ped(from) && is.pedList(to))
    return(lapply(to, .transferMarkersSimple, from = from))
  if (is.pedList(from) && is.ped(to)) {

    targetLabs = to$LABELS

    # start by transferring markers from the first in 'from'
    res = .transferMarkersSimple(from[[1]], to)
    b = as.matrix(res)

    # loop over the remaining and transfer
    for (from in from[-1]) {
      sourceLabs = from$LABELS
      shared.ids = intersect(sourceLabs, targetLabs)
      if (length(shared.ids) == 0)
        next
      a = as.matrix(from, include.attr = FALSE)
      b[match(shared.ids, targetLabs), -(1:4)] = a[match(shared.ids, sourceLabs), -(1:4)]
    }
    y = restore_ped(b)
    return(y)
  }
  if (is.pedList(from) && is.pedList(to))
    return(lapply(to, transferMarkers, from = from))
}


.transferMarkersSimple = function(from, to) {
  stopifnot(is.ped(from), is.ped(to))
  if (!hasMarkers(from)) {
    warning("No markers to transfer.")
    return(to)
  }

  # Identify shared individuals
  sourceLabs = from$LABELS
  targetLabs = to$LABELS
  shared.ids = intersect(sourceLabs, targetLabs)

  # remove prior markers in `to`
  to$markerdata = NULL

  # Prepare transfer
  a = as.matrix(from)
  b = as.matrix(to)
  b.attrs = attributes(b)

  # Transfer alleles: create empty matrix; copy rows of shared indivs
  allelematrix = a[, -(1:4), drop = FALSE]
  allelematrix.new = matrix(0L, ncol = ncol(allelematrix), nrow = pedsize(to))
  allelematrix.new[match(shared.ids, targetLabs), ] = allelematrix[match(shared.ids, sourceLabs), ]
  b = cbind(b, allelematrix.new)

  # Transfer marker attributes
  b.attrs$markerattr = attr(a, 'markerattr')

  restore_ped(b, attrs = b.attrs)
}

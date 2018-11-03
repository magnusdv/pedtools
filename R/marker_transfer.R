#' Transfer marker data
#'
#' Transfer marker data between pedigrees. Any markers attached to the target
#' are overwritten.
#'
#' @param from a `ped` or `singleton` object, or a list of such objects.
#' @param to a `ped` or `singleton` object, or a list of such objects.
#' @param ids a vector of ID labels. The indicated individuals must be present
#'   in both pedigrees. By default, genotypes are transferred between all shared
#'   individuals.
#' @param erase a logical. If `TRUE` (default), all markerdata in `to` are
#'   erased prior to transfer, and the markers in `from` are transferred
#'   including all locus annotations (frequencies a.s.o.). If `FALSE` it is
#'   assumed that the pedigree `to` already has attached markers with the same
#'   names as `from`. In this case only the genotypes of the `ids` individuals
#'   are modified; genotypes for other pedigree members are untouched, as are
#'   the marker annotations.
#' @return A `ped` object (or a list of such) similar to `to`, but where all
#'   individuals also present in `from` have marker genotypes copied over.  Any
#'   previous marker data is erased.
#'
#' @examples
#'
#' x = nuclearPed(fa = "father", mo = "mother", children = "boy")
#' m = marker(x, father = 1:2, mother = 1, boy = 1:2)
#' x = setMarkers(x, m)
#'
#' y = list(singleton("father"), nuclearPed(mo = "mother", children = "boy"))
#'
#' # By default all common individuals are transferred
#' transferMarkers(x, y)
#'
#' # Transfer data for the boy only
#' transferMarkers(x, y, ids = "boy")
#'
#' ### Transfer without erasing annotations and others genotypes
#' z = nuclearPed(children = "boy")
#' z = setMarkers(z, marker(z, '1' = c(2,2), alleles = 1:2, afreq = c(.1, .9)))
#' z2 = transferMarkers(x, z, ids = "boy", erase = FALSE)
#' z2
#' # Frequencies are not transferred
#' afreq(z2, 1)
#'
#' @export
transferMarkers = function(from, to, ids = NULL, erase = TRUE) {

  if (is.ped(from) && is.ped(to)) {
    return(.transferMarkersSimple(from, to, ids = ids, erase = erase))
  }

  if (is.ped(from) && is.pedList(to)) {

    if(is.null(ids)) # slightly cumbersome in order to catch (only relevant) repeats
      ids = unlist(lapply(to, function(p) intersect(labels(from), labels(p))))
    if((dup <- anyDuplicated(ids)) > 0)
      stop2("ID label is not unique: ", ids[dup])

    to = lapply(to, function(comp) {
      ids_comp = intersect(ids, labels(comp))
      .transferMarkersSimple(from, comp, ids = ids_comp, erase = erase)
    })
    return(to)
  }

  if (is.pedList(from) && is.ped(to)) {

    if(is.null(ids)) # slightly cumbersome in order to catch (only relevant) repeats
      ids = unlist(lapply(from, function(p) intersect(labels(p), labels(to))))

    if((dup <- anyDuplicated(ids)) > 0)
      stop2("ID label is not unique: ", ids[dup])

    # Transfer from first component
    ids1 = intersect(ids, labels(from[[1]]))
    to = .transferMarkersSimple(from[[1]], to, ids = ids1, erase = erase)

    # Transfer from remaining comps, with erase = FALSE
    for (comp in from[-1]) {
      ids_comp = intersect(ids, labels(comp))
      if (!length(ids_comp))
        next
      to = .transferMarkersSimple(comp, to, ids = ids_comp, erase = FALSE)
    }
    return(to)
  }

  if (is.pedList(from) && is.pedList(to))
    return(lapply(to, transferMarkers, from = from, ids = ids, erase = erase))
}


.transferMarkersSimple = function(from, to, ids = NULL, erase = TRUE) {
  stopifnot(is.ped(from), is.ped(to))
  M = nMarkers(from)

  if (M == 0) {
    warning("No markers to transfer.")
    return(to)
  }

  # If `ids` not given: use all shared individuals
  if(is.null(ids)) {
    ids = intersect(labels(from),  labels(to))
  }

  # Internal indices (=matrix row numbers below); also catches wrong IDs!
  idx_from = internalID(from, ids)
  idx_to = internalID(to, ids)

  # If "erase": remove all markers in `to`.
  # Otherwise: select all markers in `from`, in same order
  if(erase)
    to = setMarkers(to, NULL)
  else
    to = selectMarkers(to, markers = name(from, 1:M))

  # Prepare transfer
  a = as.matrix(from)
  b = as.matrix(to)
  b.attrs = attributes(b)

  # Allele matrices
  if(erase) {
    # Add empty allele matrix
    b = cbind(b, matrix(0L, ncol = 2*M, nrow = pedsize(to)))

    # Transfer loci annotations
    b.attrs$markerattr = attr(a, 'markerattr')
  }

  # Transfer alleles
  b[idx_to, -(1:4)] = a[idx_from, -(1:4)]

  # Restore `to` and return
  restore_ped(b, attrs = b.attrs)
}

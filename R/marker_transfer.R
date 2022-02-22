#' Transfer marker data
#'
#' Transfer marker data between pedigrees. Any markers attached to the target
#' are overwritten.
#'
#' By default, genotypes are transferred between all individuals present in both
#' pedigrees.
#'
#' @param from A `ped` or `singleton` object, or a list of such objects.
#' @param to A `ped` or `singleton` object, or a list of such objects.
#' @param ids A vector of ID labels. This should be used only if the individuals
#'   have the same name in both pedigrees; otherwise use `idsFrom` and `idsTo`
#'   instead.
#' @param idsFrom,idsTo Vectors of equal length, denoting source individuals (in
#'   the `from` pedigree) and target individuals (in the `to` pedigree),
#'   respectively.
#' @param erase A logical. If `TRUE` (default), all markers attached to `to` are
#'   erased prior to transfer, and new marker objects are created with the same
#'   attributes as in `from`. If `FALSE` no new marker objects are attached to
#'   `to`. Only the genotypes of the `ids` individuals are modified, while
#'   genotypes for other pedigree members - and marker attributes - remain
#'   untouched.
#' @param matchNames A logical, only relevant if `erase = FALSE`. If `matchNames
#'   = TRUE` (default) marker names are used to ensure genotypes are transferred
#'   into the right markers, The output contains only markers present in `from`,
#'   in the same order. (An error is raised if the markers are not named.)
#' @param checkSex A logical. If TRUE, it is checked that `fromIds` and `toIds`
#'   have the same sex. Default: FALSE.
#'
#' @return A `ped` object (or a list of such) similar to `to`, but where all
#'   individuals also present in `from` have marker genotypes copied over.  Any
#'   previous marker data is erased.
#'
#' @examples
#'
#' x = nuclearPed(fa = "A", mo = "B", child = "C")
#' x = addMarker(x, A = "1/2", B = "1/1", C = "1/2", name = "M1")
#'
#' y = list(singleton("A"), nuclearPed(fa = "D", mo = "B", child = "C"))
#'
#' # By default all common individuals are transferred
#' transferMarkers(x, y)
#'
#' # Transfer data for the boy only
#' transferMarkers(x, y, ids = "C")
#'
#' # Transfer without first erasing the target markers
#' z = nuclearPed(fa = "A", mo = "B", child = "C")
#' z = addMarker(z, A = "1/1", alleles = 1:2, name = "M1")
#'
#' transferMarkers(x, z, ids = "C", erase = FALSE)
#' transferMarkers(x, z, ids = "C", erase = TRUE) # note the difference
#'
#' @export
transferMarkers = function(from, to, ids = NULL, idsFrom = ids, idsTo = ids,
                           erase = TRUE, matchNames = TRUE, checkSex = FALSE) {

  allFrom = unlist(labels(from))
  allTo = unlist(labels(to))

  # If ids not given: use all shared individuals
  if(is.null(idsFrom) && is.null(idsTo))
    idsFrom = idsTo = intersect(allFrom, allTo)
  else if(length(idsFrom) != length(idsTo))
    stop2(sprintf("Mismatch in transfer individuals:\n `idsFrom` = %s\n `idsTo` = %s",
                  toString(idsFrom), toString(idsTo)))

  # Check for duplicates
  if(dup <- anyDuplicated(allFrom[allFrom %in% idsFrom]))
    stop2("Non-unique ID label in source ped: ", allFrom[dup])
  if(dup <- anyDuplicated(allTo[allTo %in% idsTo]))
    stop2("Non-unique ID label in target ped: ", allTo[dup])

  if(checkSex) {
    sexFrom = getSex(from, idsFrom)
    sexTo = getSex(to, idsTo)
    if(any(bad <- sexFrom > 0 & sexTo > 0 & sexFrom != sexTo)) {
      ss = c("male", "female")
      mess = sprintf(" '%s' (%s)  -->  '%s' (%s)",
                     idsFrom[bad], ss[sexFrom[bad]], idsTo[bad], ss[sexTo[bad]])
      stop2(paste0(c("Sex mismatch", mess), collapse = "\n"))
    }
  }

  if (is.ped(from) && is.ped(to)) {
    return(.transferSimple(from, to, idsFrom = idsFrom, idsTo = idsTo,
                           erase = erase, matchNames = matchNames))
  }

  if (is.ped(from) && is.pedList(to)) {

    to = lapply(to, function(comp) {
      idx = which(idsTo %in% labels(comp))
      .transferSimple(from, comp, idsFrom = idsFrom[idx], idsTo[idx],
                      erase = erase, matchNames = matchNames)
    })
    return(to)
  }

  if (is.pedList(from) && is.ped(to)) {

    # Transfer from first component
    idx1 = which(idsFrom %in% labels(from[[1]]))
    to = .transferSimple(from[[1]], to, idsFrom = idsFrom[idx1], idsTo = idsTo[idx1],
                         erase = erase, matchNames = matchNames)

    # Transfer from remaining comps, with erase = FALSE and matchNames = FALSE
    for (comp in from[-1]) {
      idx = which(idsFrom %in% labels(comp))
      to = .transferSimple(comp, to, idsFrom = idsFrom[idx], idsTo = idsTo[idx],
                           erase = FALSE, matchNames = FALSE)
    }
    return(to)
  }

  if (is.pedList(from) && is.pedList(to)) {
    to = lapply(to, function(comp) {
      idx = which(idsTo %in% labels(comp))
      transferMarkers(from, comp, idsFrom = idsFrom[idx], idsTo[idx],
                      erase = erase, matchNames = matchNames)
    })
    return(to)
  }
}

# Transfer between single ped objects
.transferSimple = function(from, to, idsFrom, idsTo, erase = TRUE, matchNames = TRUE) {
  stopifnot2(is.ped(from), is.ped(to))
  M = nMarkers(from)

  if (M == 0) {
    warning("No markers to transfer.")
    return(to)
  }

  # Internal indices (=matrix row numbers below); also catches wrong IDs!
  idx_from = internalID(from, idsFrom)
  idx_to = internalID(to, idsTo)

  # If "erase": remove all markers in `to`.
  # Otherwise: select all markers in `from`, in same order
  if(erase)
    to = setMarkers(to, NULL)
  else if(matchNames) {
    mnames = name(from, 1:M)
    if(anyNA(mnames))
      stop2("Source pedigree contains unnamed markers. If you are sure the markers match, use `matchNames = FALSE`")
    to = selectMarkers(to, markers = mnames)
  }
  else {
    if(nMarkers(to) != M)
      stop2("Argument `matchNames = FALSE` is used, but source and target have different number of markers")
  }

  # Prepare transfer
  a = as.matrix(from)
  b = as.matrix(to)
  b.attrs = attributes(b)

  # Allele matrices
  if(erase) {
    # Add empty allele matrix
    b = cbind(b, matrix(0L, ncol = 2*M, nrow = pedsize(to)))

    # Transfer locus attributes
    b.attrs$markerattr = attr(a, 'markerattr')
  }

  # Transfer alleles
  b[idx_to, -(1:4)] = a[idx_from, -(1:4)]

  # Restore `to` and return
  restorePed(b, attrs = b.attrs, validate = FALSE)
}


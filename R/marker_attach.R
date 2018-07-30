#' Attach markers to pedigrees
#'
#' In many applications it is useful to _attach_ markers to their associated
#' `ped` object. In particular for bigger projects with many markers, this makes
#' it easier to manipulate the dataset as a unit. The function `setMarkers()`
#' replaces all existing markers with the supplied ones, while `addMarkers()`
#' appends them to any existing ones.
#'
#' @param x A `ped` object
#' @param m Either a single `marker` object, a list of `marker` objects, or a
#'   data.frame or matrix.
#' @param annotations A list of marker annotations.
#'
#' @return A `ped` object.
#' @examples
#' x = singleton(1)
#' m1 = marker(x, '1'=1:2)
#' m2 = marker(x, '1'='a')
#'
#' x = setMarkers(x, m1)
#' x = addMarkers(x, m2)
#' x
#'
#' # Reversing the order of the markers
#' x = setMarkers(x, list(m2, m1))
#' x
#'
#' @name marker_attach
NULL

#' @rdname marker_attach
#' @export
setMarkers = function(x, m, annotations = NULL) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  if (is.null(m)) {
    x['markerdata'] = list(NULL)
    return(x)
  }

  mlist = NULL
  if (is.marker(m))
    mlist = list(m)
  else if (is.markerList(m))
    mlist = m
  else if (!(is.data.frame(m) || is.matrix(m)))
    stop("Argument 'm' must be either:\n",
         " * a single `marker` object`\n",
         " * a list of `marker` objects\n",
         " * a data.frame or matrix.", call.=F)

  # If markerlist, attach to x and return
  if(!is.null(mlist)) {
    class(mlist) = "markerList"
    checkConsistency(x, mlist)
    x$markerdata = mlist
    return(x)
  }

  m = as.matrix(m)
  nc = ncol(m)
  nr = nrow(m)
  if(nr != pedsize(x))
    stop("Incompatible input. Pedigree has size ", pedsize(x),
         " but allele matrix has ", nr, " rows.", calls.=FALSE)

  # If input has 1 genotype per column, split alleles to separate columns
  is_merged = is.character(m) && grepl("/", m[1,1])
  if (is_merged) {
    splitvec = unlist(strsplit(m, "/", fixed = T))
    msplit = matrix(0, ncol = 2 * nc, nrow = nr)
    msplit[, 2 * seq_len(nc) - 1] = splitvec[2 * seq_len(nc * nr) - 1]
    msplit[, 2 * seq_len(nc)] = splitvec[2 * seq_len(nc * nr)]
    m = msplit
  }

  if (ncol(m) %% 2 != 0)
    stop2("Uneven number of marker allele columns")

  nMark = ncol(m)/2

  if (!is.null(annotations)) {
    if (length(annotations) == 2 && !is.null(names(annotations)))
      annotations = rep(list(annotations), nMark)  # if given attrs for a single marker
    else if (length(annotations) != nMark)
      stop2(sprintf("Length of annotation list (%d) does not equal number of markers (%d)",
                    length(annotations), nMark))

    mlist = lapply(1:nMark, function(i) {
      if (is.null(attribs <- annotations[[i]]))
        return(NULL)
      mi = m[, c(2 * i - 1, 2 * i), drop = FALSE]
      do.call(.createMarkerObject, c(list(matr = mi), attribs))
    })
  }
  else {
    mlist = lapply(1:nMark, function(i) {
      mi = m[, c(2 * i - 1, 2 * i), drop = FALSE]
      .createMarkerObject(mi)
    })
  }
  mlist[sapply(mlist, is.null)] = NULL
  class(mlist) = "markerList"
  x$markerdata = mlist
  x
}

#' @rdname marker_attach
#' @export
addMarkers = function(x, m, annotations = NULL) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  if (is.null(m))
    return(x)

  mlist = NULL
  if (is.marker(m))
    mlist = list(m)
  else if (is.markerList(m))
    mlist = m
  else stop2("Matrix or data.frame input is not implemented yet")

  # If markerlist, append to x and return
  if(!is.null(mlist)) {
    checkConsistency(x, mlist)
    x$markerdata = c(x$markerdata, mlist)
    class(x$markerdata) = "markerList"
    return(x)
  }
}


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

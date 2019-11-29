#' Select or remove attached markers
#'
#' Functions for manipulating markers attached to `ped` objects.
#'
#' @param x A `ped` object, or a list of such
#' @param markers Either a character vector (with marker names), a numeric
#'   vector (with marker indices), a logical (of the same length as the number
#'   of markers attached to `x`), or NULL
#' @param chroms A vector of chromosome names, or NULL
#' @param fromPos A single number or NULL
#' @param toPos A single number or NULL
#'
#' @return The return values of these functions are:
#'
#'   * `selectMarkers()`: an object identical to `x`, but where only the
#'   indicated markers are kept
#'
#'   * `removeMarkers()`: an object identical to `x`, but where the indicated
#'   markers are removed
#'
#'   * `getMarkers()`: a list of `marker` objects. Note: If `x` is a list of
#'   pedigrees, the marker objects attached to the first component will be
#'   returned.
#'
#'   * `whichMarkers()`: an integer vector with indices of the indicated
#'   markers. If `x` is a list of pedigrees an error is raised unless
#'   `whichMarkers()` gives the same result for all components.
#'
#' @seealso [`setMarkers()`]
#'
#' @name marker_select
NULL

#' @rdname marker_select
#' @export
selectMarkers = function(x, markers = NULL, chroms = NULL, fromPos = NULL, toPos = NULL) {
  if(is.pedList(x)) {
    y = lapply(x, function(comp)
      selectMarkers(comp, markers = markers, chroms = chroms, fromPos = fromPos, toPos = toPos))
    return(y)
  }

  idx = whichMarkers(x, markers = markers, chroms = chroms, fromPos = fromPos, toPos = toPos)
  x$MARKERS = x$MARKERS[idx]
  x
}

#' @rdname marker_select
#' @export
getMarkers = function(x, markers = NULL, chroms = NULL, fromPos = NULL, toPos = NULL) {
  if(is.pedList(x))
    x = x[[1]]

  idx = whichMarkers(x, markers = markers, chroms = chroms, fromPos = fromPos, toPos = toPos)
  x$MARKERS[idx]
}

#' @rdname marker_select
#' @export
removeMarkers = function(x, markers = NULL, chroms = NULL, fromPos = NULL, toPos = NULL) {
  if(is.pedList(x)) {
    y = lapply(x, function(comp)
      removeMarkers(comp, markers = markers, chroms = chroms, fromPos = fromPos, toPos = toPos))
    return(y)
  }

  idx = whichMarkers(x, markers = markers, chroms = chroms, fromPos = fromPos, toPos = toPos)
  x$MARKERS[idx] = NULL
  x
}

#' @rdname marker_select
#' @export
whichMarkers = function(x, markers = NULL, chroms = NULL, fromPos = NULL, toPos = NULL) {

  if(is.pedList(x)) {
    idxList = lapply(x, function(comp)
      whichMarkers(comp, markers = markers, chroms = chroms, fromPos = fromPos, toPos = toPos))
    if(!listIdentical(idxList))
      stop2("The output of `whichMarkers()` differs between components")
    return(idxList[[1]])
  }

  # Early returns if nothing to do
  if(is.null(markers) && is.null(chroms) && is.null(fromPos) && is.null(toPos))
    return(integer(0))

  nMark = nMarkers(x)
  if(nMark == 0)
    return(integer(0))

  stopifnotSimpleVector(markers, "markers")

  if (is.logical(markers)) {
    if(length(markers) != nMark)
      stop2(sprintf("`markers` is a logical of length %d, which differs from the number of attached markers (%d)",
                    length(markers), nMark))
    idx = which(markers)
  }
  else if (is.numeric(markers)) {
    idx = as.integer(markers)

    if(any(markers != idx))
      stop2("Marker index must be integer: ", markers[markers != idx])

    outside = idx < 1 | idx > nMark
    if(any(outside))
      stop2("Marker index out of range: ", markers[outside])

  }
  else if (is.character(markers)) {
    allnames = vapply(x$MARKERS, name.marker, character(1))
    idx = match(markers, allnames)

    if(anyNA(idx))
      stop2("Unknown marker name: ", markers[is.na(idx)])
  }
  else if(is.null(markers))
    idx = seq_len(nMark)
  else
    stop2("Argument `markers` must be a numeric, a character, or NULL")

  # Return if already empty
  if (length(idx) == 0)
    return(idx)

  ### chrom
  if (!is.null(chroms)) {
    chrom_attr = vapply(x$MARKERS[idx], chrom.marker, character(1))
    idx = idx[chrom_attr %in% chroms]
  }

  ### pos
  if (!is.null(fromPos)) {
    pos_attr = vapply(x$MARKERS[idx], posMb.marker, 1)
    idx = idx[pos_attr >= fromPos]
  }
  if (!is.null(toPos)) {
    pos_attr = vapply(x$MARKERS[idx], posMb.marker, 1)
    idx = idx[pos_attr <= toPos]
  }
  idx
}


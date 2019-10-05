#' Select or remove attached markers
#'
#' Functions for manipulating markers attached to a `ped` object.
#'
#' @param x A `ped` object
#' @param markers Either a character vector (with marker names), a numeric
#'   vector (with marker indices), or NULL
#' @param chroms A vector of chromosome names, or NULL
#' @param fromPos A single number or NULL
#' @param toPos A single number or NULL
#'
#' @return The return values of these functions are:
#'
#' * `selectMarkers()` : a `ped` object where only the indicated markers are kept
#' * `removeMarkers()` : a `ped` object where the indicated markers are removed
#' * `getMarkers()` : a list of `marker` objects
#' * `whichMarkers()` : an integer vector with indices of the indicated markers
#'
#' @seealso [`setMarkers()`]
#'
#' @name marker_select
NULL

#' @rdname marker_select
#' @export
selectMarkers = function(x, markers = NULL, chroms = NULL, fromPos = NULL, toPos = NULL) {
  idx = whichMarkers(x, markers=markers, chroms=chroms, fromPos=fromPos, toPos=toPos)
  x$MARKERS = x$MARKERS[idx]
  x
}

#' @rdname marker_select
#' @export
getMarkers = function(x, markers = NULL, chroms = NULL, fromPos = NULL, toPos = NULL) {
  idx = whichMarkers(x, markers=markers, chroms=chroms, fromPos=fromPos, toPos=toPos)
  x$MARKERS[idx]
}

#' @rdname marker_select
#' @export
removeMarkers = function(x, markers = NULL, chroms = NULL, fromPos = NULL, toPos = NULL) {
  idx = whichMarkers(x, markers=markers, chroms=chroms, fromPos=fromPos, toPos=toPos)
  x$MARKERS[idx] = NULL
  x
}

#' @rdname marker_select
#' @export
whichMarkers = function(x, markers = NULL, chroms = NULL, fromPos = NULL, toPos = NULL) {

  # Early returns if nothing to do
  if(is.null(markers) && is.null(chroms) && is.null(fromPos) && is.null(toPos))
    return(integer(0))
  if(!hasMarkers(x))
    return(integer(0))

  stopifnotSimpleVector(markers, "markers")

  if (is.numeric(markers)) {
    idx = as.integer(markers)

    if(any(markers != idx))
      stop2("Marker index must be integer: ", markers[markers != idx])

    outside = idx < 1 | idx > nMarkers(x)
    if(any(outside))
      stop2("Marker index out of range: ", markers[outside])

  }
  else if (is.character(markers)) {
    allnames = vapply(x$MARKERS, name, character(1))
    idx = match(markers, allnames)

    if(anyNA(idx))
      stop2("Unknown marker name: ", markers[is.na(idx)])
  }
  else if(is.null(markers))
    idx = seq_len(nMarkers(x))
  else
    stop2("Argument `markers` must be a numeric, a character, or NULL")

  # Return if already empty
  if (length(idx) == 0)
    return(idx)

  ### chrom
  if (!is.null(chroms)) {
    chrom_attr = vapply(x$MARKERS[idx], chrom, character(1))
    idx = idx[chrom_attr %in% chroms]
  }

  ### pos
  if (!is.null(fromPos)) {
    pos_attr = vapply(x$MARKERS[idx], posMb, 1)
    idx = idx[pos_attr >= fromPos]
  }
  if (!is.null(toPos)) {
    pos_attr = vapply(x$MARKERS[idx], posMb, 1)
    idx = idx[pos_attr <= toPos]
  }
  idx
}


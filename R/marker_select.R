#' Select or remove marker objects
#'
#' Get the index (position in the markerdata list) of one or several markers
#' with specified properties.
#'
#' @param x A `ped` object
#' @param markers Either a character vector (with marker names), a numeric
#'   vector (with marker indices), or NULL.
#' @param chrom An integer vector, or NULL
#' @param fromPos A single number or NULL
#' @param toPos A single number or NULL
#'
#' @return Each of `selectMarkers()` and `removeMarkers()` returns a `ped`
#'   object, where the specified markers are kept/removed. The function
#'   `getMarkers()` returns only the `markerList` object, while `whichMarkers()`
#'   returns an integer vector with the indicies (position in `x$markerdata`) of
#'   the specified markers. NULL arguments are skipped, so `whichMarkers(x)`
#'   will return `seq_len(nMarkers(x))` (i.e. all markers).
#' @seealso [`setMarkers()`]
#'
#' @export
selectMarkers = function(x, markers = NULL, chrom = NULL, fromPos = NULL, toPos = NULL) {
  idx = whichMarkers(x, markers=markers, chrom=chrom, fromPos=fromPos, toPos=toPos)
  x$markerdata = x$markerdata[idx]
  x
}

#' @export
#' @rdname selectMarkers
getMarkers = function(x, markers = NULL, chrom = NULL, fromPos = NULL, toPos = NULL) {
  idx = whichMarkers(x, markers=markers, chrom=chrom, fromPos=fromPos, toPos=toPos)
  x$markerdata[idx]
}

#' @export
#' @rdname selectMarkers
removeMarkers = function(x, markers = NULL, chrom = NULL, fromPos = NULL, toPos = NULL) {
  idx = whichMarkers(x, markers=markers, chrom=chrom, fromPos=fromPos, toPos=toPos)
  x$markerdata[idx] = NULL
  x
}

#' @export
#' @rdname selectMarkers
whichMarkers = function(x, markers = NULL, chrom = NULL, fromPos = NULL, toPos = NULL) {

  # Early returns if nothing to do
  if(is.null(markers) && is.null(chrom) && is.null(fromPos) && is.null(toPos))
    return(integer(0))
  if(!hasMarkers(x))
    return(integer(0))

  # Argument `markers` is either numeric, character or NULL
  if (is.numeric(markers)) {
    idx = as.integer(markers)

    if(any(markers != idx)) {
      stop("Marker index must be integer: ",
           paste(markers[markers != idx], collapse=","), call.=F)
    }
    outside = idx < 1 | idx > nMarkers(x)
    if(any(outside))
      stop("Marker index out of range: ", paste(markers[outside], collapse=", "), call.=F)

  }
  else if (is.character(markers)) {
    allnames = vapply(x$markerdata, name, character(1))
    idx = match(markers, allnames)

    if(anyNA(idx)) {
      na_name = markers[is.na(idx)]
      stop("Unknown marker name: ", paste(na_name, collapse=", "), call.=F)
    }
  }
  else if(is.null(markers))
    idx = seq_len(nMarkers(x))
  else
    stop("Argument `markers` must be a numeric, a character, or NULL", call. =F)

  # Return if already empty
  if (length(idx) == 0)
    return(idx)

  ### chrom
  if (!is.null(chrom)) {
    chrom_attr = vapply(x$markerdata[idx], chrom, character(1))
    idx = idx[chrom_attr %in% chrom]
  }

  ### pos
  if (!is.null(fromPos)) {
    pos_attr = vapply(x$markerdata[idx], posMb, 1)
    idx = idx[pos_attr >= fromPos]
  }
  if (!is.null(toPos)) {
    pos_attr = vapply(x$markerdata[idx], posMb, 1)
    idx = idx[pos_attr <= toPos]
  }
  idx
}


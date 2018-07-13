#' Select or remove marker objects
#'
#' Get the index (position in the markerdata list) of one or several markers
#' with specified properties.
#'
#' @param x A `ped` object
#' @param idx A numeric index vector specifying the markers to be
#'   selected/removed.
#' @param markernames A character vector with marker names, or NULL
#' @param chroms An integer vector, or NULL
#' @param fromPos A single number or NULL
#' @param toPos A single number or NULL
#'
#' @return Each of `selectMarkers()` and `removeMarkers()` returns a `ped`
#'   object, where the specified markers are attached/removed. The function
#'   `getMarkers()` returns only the `markerList` object, while `whichMarkers()`
#'   returns an integer vector with the indicies (position in `x$markerdata`) of
#'   the specified markers. NULL arguments are skipped, so `whichMarkers(x)`
#'   will return `seq_len(nMarkers(x))` (i.e. all markers).
#' @seealso [`setMarkers()`]
#'
#' @export
selectMarkers = function(x, idx = NULL, markernames = NULL, chroms = NULL, fromPos = NULL, toPos = NULL) {
  if(is.null(idx)) {
    idx = whichMarkers(x, markernames = markernames, chroms = chroms, fromPos = fromPos, toPos = toPos)
  }
  assert_that(is.numeric(idx), all(idx == as.integer(idx)), all(idx <= nMarkers(x))) # TRUE also when empty
  if (length(idx) == 0)
    x$markerdata = NULL
  else
    x$markerdata = x$markerdata[idx]
  x
}


#' @export
#' @rdname selectMarkers
getMarkers = function(x, idx = NULL, markernames = NULL, chroms = NULL, fromPos = NULL, toPos = NULL) {
  if(is.null(idx)) {
    idx = whichMarkers(x, markernames = markernames, chroms = chroms, fromPos = fromPos, toPos = toPos)
  }
  assert_that(is.numeric(idx), all(idx == as.integer(idx)), all(idx <= nMarkers(x))) # TRUE also when empty
  x$markerdata[idx]
}

#' @export
#' @rdname selectMarkers
removeMarkers = function(x, idx = NULL, markernames = NULL, chroms = NULL, fromPos = NULL, toPos = NULL) {
  if (is.null(markers))
    markers = whichMarkers(x, markernames, chroms, fromPos, toPos)
  if (is.null(markers) || length(markers) == 0)
    return(x)
  m = x$markerdata
  m[markers] = NULL
  setMarkers(x, m)
}

#' @export
#' @rdname selectMarkers
whichMarkers = function(x, markernames = NULL, chroms = NULL, fromPos = NULL, toPos = NULL) {
  mnos = seq_len(nMarkers(x))
  if (!is.null(markernames)) {
    if (length(markernames) == 0)
      return(numeric(0))

    name_attrs = unlist(lapply(x$markerdata, function(m) attr(m, "name")))
    mnos = mnos[match(markernames, name_attrs, nomatch = 0)]
  }
  if (!is.null(chroms)) {
    if (length(chroms) == 0)
      return(numeric(0))

    chrom_attrs = unlist(lapply(x$markerdata[mnos], function(m) attr(m, "chrom")))
    mnos = mnos[chrom_attrs %in% chroms]
  }
  if (!is.null(fromPos))
    mnos = mnos[unlist(lapply(x$markerdata[mnos], function(m) attr(m, "pos"))) >= fromPos]
  if (!is.null(toPos))
    mnos = mnos[unlist(lapply(x$markerdata[mnos], function(m) attr(m, "pos"))) <= toPos]
  mnos
}


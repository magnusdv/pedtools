

# Check if a mutation matrix - or a list of such - is diagonal (= trivial)
trivialMut = function(mut) {
  if(is.list(mut))
    return(all(vapply(mut, trivialMut, logical(1))))

  all(diag(mut) == 1)
}

#' Test if something is a marker
#'
#' Functions for testing if something is a `marker` object, or a list of such objects.
#'
#' @param x Any object
#'
#' @return A logical
#' @export
is.marker = function(x) {
  inherits(x, "marker")
}

#' @export
#' @rdname is.marker
is.markerList = function(x) {
  inherits(x, "markerList") || (is.list(x) && all(sapply(x, is.marker)))
}


#' The number of markers attached to a pedigree
#'
#' @param x A `ped` object or a list of such (se Value).
#' @return The function `nMarkers` returns the number of marker objects attached
#'   to `x`. If `x` is a list of pedigrees, an error is raised unless all of
#'   them have the same number of markers.
#'
#'   The function `hasMarkers` returns TRUE if `nMarkers(x) > 0`.
#'
#' @export
nMarkers = function(x) {
  if(is.ped(x))
    return(length(x$MARKERS))
  else if(is.pedList(x)) {
    nvec = vapply(x, function(comp) length(comp$MARKERS), 1L)
    if(!listIdentical(nvec))
      stop2("The pedigree components have different number of markers attached")
    return(nvec[[1]])
  }
  stop2("Input to `nMarkers()` must be a `ped` object or a list of such")
}

#' @export
#' @rdname nMarkers
hasMarkers = function(x) {
  nMarkers(x) > 0
}


checkConsistency = function(x, mlist) {
  wrongSize = unlist(lapply(mlist, nrow) != pedsize(x))
  if(any(wrongSize)) {
    erri = which(wrongSize)[1]
    errsize = nrow(mlist[[erri]])
    stop("Incompatible input: Pedigree has size ", pedsize(x),
         " but marker ", erri, " has ", errsize, " rows", call. = FALSE)
  }
  #TODO: check loop breakers, ped labels, sex
  return(TRUE)
}



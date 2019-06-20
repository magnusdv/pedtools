
#' Check if a marker allows mutations
#'
#' @param marker A `marker` object
#'
#' @return Returns TRUE if the `mutmod` attribute of the input is non-NULL and differs from the identity matrix.
#' @export
allowsMutations = function(marker) {
  mut = mutmod(marker)
  !is.null(mut) && !trivialMut(mut)
}

# Check if a mutation matrix - or a list of such - is diagonal (= trivial)
trivialMut = function(mut) {
  if(is.list(mut))
    return(all(vapply(mut, trivialMut, logical(1))))
  all(diag(mut) == 1)
}

#' Number of marker alleles
#'
#' @param m An object of class [marker].
#'
#' @return A positive integer.
#' @export
#'
#' @examples
#' x = nuclearPed(1)
#' m = marker(x)
#' nAlleles(m)
#'
nAlleles = function(m) {
  if(!is.marker(m)) stop2("Input is not a `marker` object")
  length(attr(m, 'alleles'))
}


#' Test if something is a marker or a markerList
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

#' Test if a marker is on the X chromosome
#'
#' Tests if the `chrom` attribute of a marker is 23.
#'
#' @param x A marker object. (If `x` is not a marker object, the function
#'   returns FALSE.)
#' @return TRUE or FALSE.
#' @export
isXmarker = function(x) {
  chr = chrom(x)
  isTRUE(chr == "X" || chr == "23")
}

#' The number of markers attached to a pedigree
#'
#' @param x A `ped` object.
#' @return The function `nMarkers` returns the number of marker objects attached
#'   to `x`. The function `hasMarkers` returns TRUE if `nMarkers(x) > 0`.
#'
#' @export
nMarkers = function(x) {
  length(x$markerdata)
}

#' @export
#' @rdname nMarkers
hasMarkers = function(x) {
  !is.null(x$markerdata) && nMarkers(x) > 0
}


checkConsistency = function(x, mlist) {
  wrongSize = unlist(lapply(mlist, nrow) != pedsize(x))
  if(any(wrongSize)) {
    erri = which(wrongSize)[1]
    errsize = nrow(mlist[[erri]])
    stop("Incompatible input: Pedigree has size ", pedsize(x),
         " but marker ", erri, " has ", errsize, " rows", call.=FALSE)
  }
  #TODO: check loop breakers, ped labels, sex
  return(TRUE)
}



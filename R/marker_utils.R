
#' Check if a marker allows mutatations
#'
#' @param marker A `marker` object
#'
#' @return Returns TRUE if the `mutmat` attribute of marker is non-NULL and differs from the identity matrix.
#' @export
allowsMutations = function(marker) {
  mutmat = attr(marker, 'mutmat')
  if(is.null(mutmat))
    return(FALSE)
  n = nAlleles(marker)
  male = mutmat$male
  female = mutmat$female
  if(all(diag(male) == 1) && all(diag(female) == 1))
    return(FALSE)
  return(TRUE)
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
  assert_that(is.marker(m))
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

assertthat::on_failure(is.marker) = function (call, env) {
  paste0(deparse(call$x), " is not a `marker` object")
}

#' @export
#' @rdname is.marker
is.markerList = function(x) {
  inherits(x, "markerList") || (is.list(x) && all(sapply(x, is.marker)))
}

assertthat::on_failure(is.markerList) = function (call, env) {
  paste0(deparse(call$x), " is not a list of `marker` objects")
}

#' Test if a marker is on the X chromosome
#'
#' Tests if the `chrom` attribute of a marker is 23.
#'
#' @param x A marker object. (If `x` is not a marker object, the function
#'   returns FALSE.)
#' @return TRUE or FALSE.
#' @export
is_Xmarker = function(x) {
  assert_that(is.marker(x))
  isTRUE(attr(x, 'chrom') == 23)
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
  wrongSize = unlist(lapply(mlist, nrow) != pedSize(x))
  if(any(wrongSize)) {
    erri = which(wrongSize)[1]
    errsize = nrow(mlist[[erri]])
    stop("Incompatible input: Pedigree has size ", pedSize(x),
         " but marker ", erri, " has ", errsize, " rows", call.=FALSE)
  }
  #TODO: check loop breakers
  return(TRUE)
}



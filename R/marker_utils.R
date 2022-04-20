

# Check if a mutation matrix - or a list of such - is diagonal (= trivial)
trivialMut = function(mut) {
  if(is.list(mut))
    return(all(vapply(mut, trivialMut, logical(1))))

  all(diag(mut) == 1)
}

#' Test if something is a marker
#'
#' Functions for testing if something is a `marker` object, or a list of such
#' objects.
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
#' @param x A `ped` object or a list of such (see Value).
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

checkDupNames = function(x) {
  if(is.ped(x))
    mlist = x$MARKERS
  else if(is.pedList(x))
    mlist = x[[1]]$MARKERS
  else
    mlist = x

  mnames = unlist(lapply(mlist, attr, "name"))
  if(dup <- anyDuplicated.default(mnames, incomparables = NA))
    stop2("Duplicated marker name: ", mnames[dup])
}

checkConsistency = function(x, mlist) {
  checkDupNames(mlist)

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


#' Add allele
#'
#' Extends the allele set of a marker attached to a pedigree, by adding a single
#' allele.
#'
#' @param x A `ped` object or a list of such.
#' @param marker The name or index of a marker attached to `x`.
#' @param allele The name of the new allele.
#' @param freq The frequency of the new allele, by default 0.001.
#' @param adjust Either "previous" or "all", indicating how the frequencies
#'   should be adjusted so that they sum to 1. If "previous" (default), the
#'   frequencies of the original alleles are multiplied with `1 - freq`. If
#'   "all", scaling is performed after adding the new allele, i.e., dividing all
#'   frequencies by `1 + freq`.
#'
#' @return A copy of `x` with modified marker attributes.
#'
#' @examples
#' x = nuclearPed() |>
#'   addMarker(geno = c(NA, NA, "b/c"), afreq = c(b = 0.5, c = 0.5))
#'
#' y = addAllele(x, marker = 1, allele = "a")
#' afreq(y, 1)
#'
#' z = addAllele(y, marker = 1, allele = "d", freq = 0.1, adjust = "all")
#' afreq(z, 1)
#'
#' @export
addAllele = function(x, marker, allele, freq = 0.001, adjust = c("previous", "all")) {

  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Allele modifications must be done for a single marker at a time")
  if(missing(allele) || length(allele) == 0)
    stop2("Argument `allele` cannot be empty")
  if(length(allele) > 1)
    stop2("Please add one allele at a time")

  old = afreq(x, marker)
  if(allele %in% names(old))
    stop2("This allele is already present")

  names(freq) = allele
  switch(match.arg(adjust),
    all = {
      newfr = c(old, freq)/(1 + freq)
    },
    previous = {
      newfr = c(old * (1 - freq), freq)
    }
  )

  setAfreq(x, marker, afreq = newfr, strict = FALSE)
}

#' Marker properties
#'
#' These functions are used to retrieve various properties of marker objects.
#' Each function accepts as input either a single `marker` object, a `ped`
#' object, or a list of `ped` objects.
#'
#' `emptyMarker()` returns TRUE for markers with no genotypes. If the input is a
#' list of pedigrees, all must be empty for the result to be TRUE.
#'
#' `nTyped()` returns the number of typed individuals for each marker. Note that
#' if the input is a list of pedigrees, the function returns the sum over all
#' components.
#'
#' `nAlleles()` returns the number of alleles of each marker.
#'
#' `isXmarker()` returns TRUE for markers whose `chrom` attribute is either "X"
#' or 23.
#'
#' `allowsMutations` returns TRUE for markers whose `mutmod` attribute is
#' non-NULL and differs from the identity matrix.
#'
#' @param x A single `marker` object or a `ped` object (or a list of such)
#' @param markers A vector of names or indices of markers attached to `x`, in
#'   the case that `x` is a `ped` object or a list of such. By default all
#'   attached markers are selected.
#' @param ... Not used.
#'
#' @return
#'
#' If `x` is a single `marker` object, the output is a vector of length 1.
#'
#' If `x` is a `ped` object, or a list of such, the output is a vector of the
#' same length as `markers` (which includes all attached markers by default),
#' reporting the property of each marker.
#'
#' @examples
#' cmp1 = nuclearPed(1)
#' cmp2 = singleton(10)
#' loc = list(alleles = 1:2)
#' x = setMarkers(list(cmp1, cmp2), locus = rep(list(loc), 3))
#'
#' #-------- nAlleles() ------------
#' # All markers have 2 alleles
#' stopifnot(identical(nAlleles(x), c(2L,2L,2L)))
#'
#' #-------- emptyMarkers() ------------
#' # Add genotype for indiv 1 at marker 1
#' genotype(x[[1]], 1, 1) = "1/2"
#'
#' # Check that markers 2 and 3 are empty
#' stopifnot(identical(emptyMarker(x), c(FALSE,TRUE,TRUE)),
#'           identical(emptyMarker(x[[1]]), c(FALSE,TRUE,TRUE)),
#'           identical(emptyMarker(x[[2]]), c(TRUE,TRUE,TRUE)),
#'           identical(emptyMarker(x, markers = c(3,1)), c(TRUE,FALSE)))
#'
#' #-------- nTyped() ------------
#' stopifnot(identical(nTyped(x), c(1L,0L,0L)))
#'
#' # Add genotypes for third marker
#' genotype(x[[1]], marker = 3, id = 1:3) = "1/1"
#' genotype(x[[2]], marker = 3, id = 10) = "2/2"
#'
#' # nTyped() returns total over all components
#' stopifnot(identical(nTyped(x), c(1L,0L,4L)))
#'
#' #-------- allowsMutations() ------------
#' # Marker 2 allows mutations
#' mutmod(x, 2) = list("prop", rate = 0.1)
#'
#' stopifnot(identical(allowsMutations(x), c(FALSE,TRUE,FALSE)),
#'           identical(allowsMutations(x, markers = 2:3), c(TRUE,FALSE)))
#'
#' #-------- isXmarker() ------------
#' # Make marker 3 X-linked
#' chrom(x[[1]], 3) = "X"
#' chrom(x[[2]], 3) = "X"
#'
#' stopifnot(identical(isXmarker(x), c(FALSE,FALSE,TRUE)))
#'
#' @name marker_prop
NULL

#' @rdname marker_prop
#' @export
emptyMarker = function(x, ...) {
  UseMethod("emptyMarker")
}

#' @rdname marker_prop
#' @export
emptyMarker.default = function(x, ...) {
  stop2("No method defined for objects of class ", class(x)[1])
}

#' @rdname marker_prop
#' @export
emptyMarker.marker = function(x, ...) {
  all(x == 0)
}

#' @rdname marker_prop
#' @export
emptyMarker.ped = function(x, markers = seq_len(nMarkers(x)), ...) {
  vapply(getMarkers(x, markers), emptyMarker.marker, logical(1))
}

#' @rdname marker_prop
#' @export
emptyMarker.list = function(x, markers = seq_len(nMarkers(x)), ...) {
  if(length(x) == 0)
    return(logical(0))
  comp_wise = lapply(x, emptyMarker.ped, markers = markers)
  Reduce(`&`, comp_wise)
}


#' @rdname marker_prop
#' @export
nTyped = function(x, ...) {
  UseMethod("nTyped")
}

#' @rdname marker_prop
#' @export
nTyped.default = function(x, ...) {
  stop2("No method defined for objects of class ", class(x)[1])
}

#' @rdname marker_prop
#' @export
nTyped.marker = function(x, ...) {
  sum(x[, 1] > 0 | x[, 2] > 0)
}

#' @rdname marker_prop
#' @export
nTyped.ped = function(x, markers = seq_len(nMarkers(x)), ...) {
  vapply(getMarkers(x, markers), nTyped.marker, integer(1))
}

#' @rdname marker_prop
#' @export
nTyped.list = function(x, markers = seq_len(nMarkers(x)), ...) {
  if(length(x) == 0)
    return(logical(0))
  comp_wise = lapply(x, nTyped.ped, markers = markers)
  Reduce(`+`, comp_wise)
}


#' @rdname marker_prop
#' @export
nAlleles = function(x, ...) {
  UseMethod("nAlleles")
}

#' @rdname marker_prop
#' @export
nAlleles.default = function(x, ...) {
  stop2("No method defined for objects of class ", class(x)[1])
}

#' @rdname marker_prop
#' @export
nAlleles.marker = function(x, ...) {
  length(attr(x, 'alleles'))
}

#' @rdname marker_prop
#' @export
nAlleles.ped = function(x, markers = seq_len(nMarkers(x)), ...) {
  vapply(getMarkers(x, markers), nAlleles.marker, integer(1))
}

#' @rdname marker_prop
#' @export
nAlleles.list = function(x, markers = seq_len(nMarkers(x)), ...) {
  if(length(x) == 0)
    return(integer(0))
  comp_wise = lapply(x, nAlleles.ped, markers = markers)
  if(!listIdentical(comp_wise))
    stop2("The output of `nAlleles()` differs between components")
  comp_wise[[1]]
}


#' @rdname marker_prop
#' @export
isXmarker = function(x, ...) {
  UseMethod("isXmarker")
}

#' @rdname marker_prop
#' @export
isXmarker.default = function(x, ...) {
  if(is.null(x)) return(FALSE)
  stop2("No method defined for objects of class ", class(x)[1])
}

#' @rdname marker_prop
#' @export
isXmarker.marker = function(x, ...) {
  chr = chrom(x)
  isTRUE(chr == "X" || chr == "23")
}

#' @rdname marker_prop
#' @export
isXmarker.ped = function(x, markers = seq_len(nMarkers(x)), ...) {
  vapply(getMarkers(x, markers), isXmarker.marker, logical(1))
}

#' @rdname marker_prop
#' @export
isXmarker.list = function(x, markers = seq_len(nMarkers(x)), ...) {
  if(length(x) == 0)
    return(logical(0))
  comp_wise = lapply(x, isXmarker.ped, markers = markers)
  if(!listIdentical(comp_wise))
    stop2("The output of `isXmarker()` differs between components")
  comp_wise[[1]]
}


#' @rdname marker_prop
#' @export
allowsMutations = function(x, ...) {
  UseMethod("allowsMutations")
}

#' @rdname marker_prop
#' @export
allowsMutations.default = function(x, ...) {
  stop2("No method defined for objects of class ", class(x)[1])
}

#' @rdname marker_prop
#' @export
allowsMutations.marker = function(x, ...) {
  mut = mutmod(x)
  !is.null(mut) || !trivialMut(mut)
}

#' @rdname marker_prop
#' @export
allowsMutations.ped = function(x, markers = seq_len(nMarkers(x)), ...) {
  vapply(getMarkers(x, markers), allowsMutations.marker, logical(1))
}

#' @rdname marker_prop
#' @export
allowsMutations.list = function(x, markers = seq_len(nMarkers(x)), ...) {
  comp_wise = lapply(x, allowsMutations.ped, markers = markers)
  if(!listIdentical(comp_wise))
    stop2("The output of `allowsMutations()` differs between components")
  comp_wise[[1]]
}

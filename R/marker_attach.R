#' Attach markers to pedigrees
#'
#' In many applications it is useful to _attach_ markers to their associated
#' `ped` object. In particular for bigger projects with many markers, this makes
#' it easier to manipulate the dataset as a unit. The function `setMarkers()`
#' replaces all existing markers with the supplied ones, while `addMarkers()`
#' appends the supplied markers to any existing ones.
#'
#' @param x A `ped` object
#' @param m Either a single `marker` object or a list of `marker` objects
#' @param allele_matrix A matrix with `pedsize(x)` rows, containing the observed
#'   alleles for one or several markers. The matrix must have either 1 or 2
#'   columns per marker. If the former, then the argument `allele_sep` must be
#'   non-NULL, and will be used to split all entries.
#' @param locus_annotations A list of lists, with annotations for each marker.
#'   See [marker()] for possible entries.
#' @param missing A single character (or coercible to one) indicating the symbol
#'   for missing alleles.
#' @param allele_sep If this is a single string, each entry of `allele_matrix`
#'   is interpreted as a genotype, and will be split by calling `str_split(...,
#'   split = allele_sep, fixed = T)`. For example, if the entries are formatted
#'   as "A/B", put `allele_sep="/"`. Default: NULL.
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
setMarkers = function(x, m = NULL, allele_matrix = NULL, locus_annotations = NULL, missing=0, allele_sep=NULL) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(is.null(m) && is.null(allele_matrix)) {
    x['markerdata'] = list(NULL)
    return(x)
  }

  if (is.marker(m))
    mlist = list(m)
  else if (is.markerList(m))
    mlist = m
  else if (is.null(m))
    mlist = allelematrix2markerlist(x, allele_matrix, locus_annotations, missing, allele_sep)
  else
    stop2("Argument `m` must be either a single `marker` object, a list of such, or NULL")

  class(mlist) = "markerList"
  checkConsistency(x, mlist)
  x$markerdata = mlist
  x
}

#' @rdname marker_attach
#' @export
addMarkers = function(x, m = NULL, allele_matrix = NULL, locus_annotations = NULL, missing=0, allele_sep=NULL) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  if (is.marker(m))
    mlist = list(m)
  else if (is.markerList(m))
    mlist = m
  else if (is.null(m))
    mlist = allelematrix2markerlist(x, allele_matrix, locus_annotations, missing, allele_sep)
  else
    stop2("Argument `m` must be either a single `marker` object, a list of such, or NULL")

  # Append to x and return
  checkConsistency(x, mlist)
  mlist = c(x$markerdata, mlist)
  class(mlist) = "markerList"
  x$markerdata = mlist
  x
}

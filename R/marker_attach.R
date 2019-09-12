#' Attach markers to pedigrees
#'
#' In many applications it is useful to _attach_ markers to their associated
#' `ped` object. In particular for bigger projects with many markers, this makes
#' it easier to manipulate the dataset as a unit. The function `setMarkers()`
#' replaces all existing markers with the supplied ones, while `addMarkers()`
#' appends the supplied markers to any existing ones.
#'
#' The most general format of `locusAttributes` a list of lists, one for each
#' marker, where possible entries in the inner lists are as follows (default
#' values in parenthesis):
#'
#' * 'alleles' : a character vector with allele labels
#'
#' * 'afreq' :  a numeric vector with allele frequencies (`rep.int(1/L, L)`,
#' where `L = length(alleles)`)
#'
#' * 'chrom' : chromosome number (NA)
#'
#' * 'posMb' : physical location in megabases (NA)
#'
#' * 'posCm' : position in centiMorgan (NA)
#'
#' * 'name' : marker name (NA)
#'
#' * 'mutmod' : mutation model, or model name (NULL)
#'
#' * 'rate' : mutation model parameter (NULL)
#'
#' If `locusAttributes` is just a single list of attributes (not a list of
#' lists), then it is repeated to match the number of markers.
#'
#' Two alternative format of `locusAttributes` are allowed: If a data.frame or
#' matrix is given, an attempt is made to interprete it as a frequency database
#' in `allelic ladder` format. Such an interpretation is also attempted if
#' `locusAttributes` is a list of named frequency vectors (where the names are
#' the allele labels). See the Examples section for an example of this using the
#' `Familias::NorwegianFrequencies` database.
#'
#' @param x A `ped` object
#' @param m Either a single `marker` object or a list of `marker` objects
#' @param alleleMatrix A matrix with `pedsize(x)` rows, containing the observed
#'   alleles for one or several markers. The matrix must have either 1 or 2
#'   columns per marker. If the former, then a `sep` string must be a given, and
#'   will be used to split all entries.
#' @param locusAttributes A list of lists, with attributes for each marker. See
#'   Details for possible attributes.
#' @param missing A single character (or coercible to one) indicating the symbol
#'   for missing alleles.
#' @param sep If this is a single string, each entry of `alleleMatrix` is
#'   interpreted as a genotype, and will be split by calling `strsplit(...,
#'   split = sep, fixed = T)`. For example, if the entries are formatted as
#'   "A/B", put `sep="/"`. Default: NULL.
#' @param ... Further arguments
#'
#' @return A `ped` object.
#' @examples
#' x = singleton(1)
#' m1 = marker(x, '1' = 1:2)
#' m2 = marker(x, '1' = 'a')
#'
#' x = setMarkers(x, m1)
#' x = addMarkers(x, m2)
#' x
#'
#' # Reversing the order of the markers
#' x = setMarkers(x, list(m2, m1))
#' x
#'
#' \dontrun{
#' # Attach empty markers corresponding to a realistic forensic database
#' require(Familias)
#' data(NorwegianFrequencies)
#' x = nuclearPed(1)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies)
#' }
#'
#' @name marker_attach
NULL

#' @rdname marker_attach
#' @export
setMarkers = function(x, m = NULL, alleleMatrix = NULL, locusAttributes = NULL, missing = 0,
                      sep = NULL, ...) {

  # Check for deprecated arguments (hackish)
  dots = list(...)
  for(arg in names(dots)) {
    if(!is.na(pmatch(arg, "allele_matrix"))) {
      warning("Argument `allele_matrix` is deprecated; use `alleleMatrix` instead")
      alleleMatrix = dots[[arg]]
    }
    if(!is.na(pmatch(arg, "locus_annotations"))) {
      warning("Argument `locus_annotations` is deprecated; use `locusAttributes` instead")
      locusAttributes = dots[[arg]]
    }
    if(!is.na(pmatch(arg, "allele_sep"))) {
      warning("Argument `allele_sep` is deprecated; use `sep` instead")
      sep = dots[[arg]]
    }
  }

  if(!is.ped(x)) stop2("Input is not a `ped` object")

  # If no data, remove all markers and return
  if(is.null(m) && is.null(alleleMatrix) && length(locusAttributes) == 0) {
    x['markerdata'] = list(NULL)
    return(x)
  }

  locusAttributes = checkLocusAttribs(locusAttributes)
  L = length(locusAttributes)

  # If only locus attributes are given, create allele matrix with all 0's.
  if(is.null(m) && is.null(alleleMatrix) && L > 0) {
    alleleMatrix = matrix(missing, nrow = pedsize(x), ncol = 2*L)
  }

  if (is.marker(m))
    mlist = list(m)
  else if (is.markerList(m))
    mlist = m
  else if (is.null(m))
    mlist = allelematrix2markerlist(x, alleleMatrix, locusAttributes, missing, sep)
  else
    stop2("Argument `m` must be either a single `marker` object, a list of such, or NULL")

  class(mlist) = "markerList"
  checkConsistency(x, mlist)
  x$markerdata = mlist
  x
}

#' @rdname marker_attach
#' @export
addMarkers = function(x, m = NULL, alleleMatrix = NULL, locusAttributes = NULL, missing = 0,
                      sep = NULL, ...) {

  # Check for deprecated arguments (hackish)
  dots = list(...)
  for(arg in names(dots)) {
    if(!is.na(pmatch(arg, "allele_matrix"))) {
      warning("Argument `allele_matrix` is deprecated; use `alleleMatrix` instead")
      alleleMatrix = dots[[arg]]
    }
    if(!is.na(pmatch(arg, "locus_annotations"))) {
      warning("Argument `locus_annotations` is deprecated; use `locusAttributes` instead")
      locusAttributes = dots[[arg]]
    }
    if(!is.na(pmatch(arg, "allele_sep"))) {
      warning("Argument `allele_sep` is deprecated; use `sep` instead")
      sep = dots[[arg]]
    }
  }

  if(!is.ped(x)) stop2("Input is not a `ped` object")

  # If no data, do nothing
  if(is.null(m) && is.null(alleleMatrix) && length(locusAttributes) == 0) {
    return(x)
  }

  locusAttributes = checkLocusAttribs(locusAttributes)
  L = length(locusAttributes)

  # If only locus attributes are given, create allele matrix with all 0's.
  if(is.null(m) && is.null(alleleMatrix) && L > 0) {
    alleleMatrix = matrix(missing, nrow = pedsize(x), ncol = 2*L)
  }

  if (is.marker(m))
    mlist = list(m)
  else if (is.markerList(m))
    mlist = m
  else if (is.null(m))
    mlist = allelematrix2markerlist(x, alleleMatrix, locusAttributes, missing, sep)
  else
    stop2("Argument `m` must be either a single `marker` object, a list of such, or NULL")

  # Append to x and return
  checkConsistency(x, mlist)
  mlist = c(x$markerdata, mlist)
  class(mlist) = "markerList"
  x$markerdata = mlist
  x
}


checkLocusAttribs = function(a) {
  if(length(a) == 0) return(a)

  attribNames = c("alleles", "afreq", "name" ,"chrom" ,"posMb", "posCm", "mutmod", "rate")

  # Format 1: List of lists
  if(is.list(a) && all(sapply(a, is.list))) {
    for(i in seq_along(a)) {
      nms = names(a[[i]])
      if(is.null(nms))
        stop2("Entry ", i, " of `locusAttributes` has no names")
      if(!all(nms %in% attribNames))
        stop2("Entry ", i, " of `locusAttributes` has illegal entries: ", setdiff(nms, attribNames))
    }
    return(a)
  }

  # Format 2: Single list of attributes
  if(is.list(a) && !is.list(a[[1]]) && !is.null(names(a)) && all(names(a) %in% attribNames)) {
    return(list(a))
  }

  # Format 3: Allelic ladder as data.frame or matrix
  if(is.data.frame(a) || is.matrix(a)) {
    return(freqDb2attribList(a, format = "allelicLadder"))
  }

  # Format 4: Frequency database as a list of frequency vectors
  if(is.list(a) && all(sapply(a, function(aa) !is.null(names(aa))))) {
    return(freqDb2attribList(a, format = "list"))
  }

}

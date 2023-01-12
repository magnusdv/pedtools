#' Attach markers to pedigrees
#'
#' In many applications it is useful to _attach_ markers to their associated
#' `ped` object. In particular for bigger projects with many markers, this makes
#' it easier to manipulate the dataset as a unit. The function `setMarkers()`
#' replaces all existing markers with the supplied ones, while `addMarkers()`
#' appends the supplied markers to any existing ones. Note that there is also
#' the function [addMarker()], which creates and attaches a single marker in one
#' go.
#'
#' The most general format of `locusAttributes` a list of lists, one for each
#' marker, where possible entries in the inner lists are as follows (default
#' values in parenthesis):
#'
#' * `alleles` : a character vector with allele labels
#'
#' * `afreq` :  a numeric vector with allele frequencies (`rep.int(1/L, L)`,
#' where `L = length(alleles)`)
#'
#' * `chrom` : chromosome number (NA)
#'
#' * `posMb` : physical location in megabases (NA)
#'
#' * `name` : marker name (NA)
#'
#' * `mutmod` : mutation model, or model name (NULL)
#'
#' * `rate` : mutation model parameter (NULL)
#'
#' If `locusAttributes` is a single list of attributes (not a list of lists),
#' then it is repeated to match the number of markers.
#'
#' #### Alternative formats of `locusAttributes`:
#'
#' * data frame or matrix. In this case an attempt is made to interpret it as a
#' frequency database in `allelic ladder` format.
#'
#' * A list of frequency vectors. All vectors should sum to 1, and be named
#' (with allele labels)
#'
#' * Shortcut for simple SNP data: The argument `locusAttributes = "snp-AB"`
#' sets all markers to be equifrequent SNPs with alleles A and B. The letters A
#' and B may be replaced by other single-character letters or numbers.
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
#'   split = sep, fixed = TRUE)`. If `alleleMatrix` contains entries with "/",
#'   this will be taken as separator by default. (To override this behaviour,
#'   put `sep = FALSE`.)
#' @param checkCons A logical. If TRUE (default), each marker is checked for
#'   consistency with `x`.
#'
#' @return A `ped` object.
#'
#' @seealso [addMarker()]
#'
#' @examples
#' x = singleton(1)
#' m1 = marker(x, `1` = "1/2")
#' m2 = marker(x, `1` = "a/b")
#'
#' # Attach to x
#' x1 = setMarkers(x, list(m1, m2))
#'
#' # Reversing the order of the markers
#' setMarkers(x, list(m2, m1))
#'
#' # Alternative syntax, adding one marker at a time
#' x2 = x |>
#'   addMarker(`1` = "1/2") |>
#'   addMarker(`1` = "a/b")
#'
#' stopifnot(identical(x1, x2))
#'
#' @name marker_attach
NULL

#' @rdname marker_attach
#' @export
setMarkers = function(x, m = NULL, alleleMatrix = NULL, locusAttributes = NULL, missing = 0,
                      sep = NULL, checkCons = TRUE) {

  # If `sep` is not given, but AM contains entries with "/", use this
  if(!is.null(alleleMatrix) && is.null(sep) && any(grepl("/", alleleMatrix, fixed = TRUE)))
    sep = "/"
  if(isFALSE(sep))
    sep = NULL

  # If pedlist input, recurse over components
  if(is.pedList(x)) {
    if(!is.null(m))
      stop2("When `x` is a list of pedigrees, argument `m` must be NULL")
    y = lapply(x, function(comp)
      setMarkers(comp, alleleMatrix = alleleMatrix, locusAttributes = locusAttributes,
                 missing = missing, sep = sep, checkCons = checkCons))
    return(y)
  }

  if(!is.ped(x))
    stop2("First argument must be a `ped` object or a list of such")

  # If no data, remove all markers and return
  if(is.null(m) && is.null(alleleMatrix) && length(locusAttributes) == 0) {
    x['MARKERS'] = list(NULL)
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

  # Check consistency with `x`
  if(checkCons)
    checkConsistency(x, mlist)

  class(mlist) = "markerList"
  x$MARKERS = mlist
  x
}


#' @rdname marker_attach
#' @export
addMarkers = function(x, m = NULL, alleleMatrix = NULL, locusAttributes = NULL, missing = 0,
                      sep = NULL, checkCons = TRUE) {

  # If `sep` is not given, but AM contains entries with "/", use this
  if(!is.null(alleleMatrix) && is.null(sep) && any(grepl("/", alleleMatrix, fixed = TRUE)))
    sep = "/"
  if(isFALSE(sep))
    sep = NULL

  # If pedlist input, recurse over components
  if(is.pedList(x)) {
    if(!is.null(m))
      stop2("When `x` is a list of pedigrees, argument `m` must be NULL")
    y = lapply(x, function(comp)
      addMarkers(comp, alleleMatrix = alleleMatrix, locusAttributes = locusAttributes,
                 missing = missing, sep = sep, checkCons = checkCons))
    return(y)
  }

  if(!is.ped(x))
    stop2("First argument must be a `ped` object or a list of such")


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

  # Check consistency with `x`
  if(checkCons)
    checkConsistency(x, mlist)

  # Append to x and return
  mlist = c(x$MARKERS, mlist)
  class(mlist) = "markerList"
  x$MARKERS = mlist
  x
}


checkLocusAttribs = function(a) {
  if(length(a) == 0)
    return(a)

  attribNames = c("alleles", "afreq", "name" ,"chrom" ,"posMb", "mutmod", "rate")

  # Shortcut for SNPs
  if(length(a) == 1 && is.character(a) && isTRUE(startsWith(tolower(a), "snp"))) {
    nch = nchar(a)
    if(!nch %in% 5:6)
      stop2("Shortcut code for SNP markers must be of the form 'snpAB' or 'snp-AB': ", a)
    a = list(alleles = strsplit(a, "")[[1]][c(nch - 1, nch)])
  }

  # Format 1: List of lists
  if(is.list(a) && all(sapply(a, is.list))) {
    for(i in seq_along(a)) {
      nms = names(a[[i]])
      if(is.null(nms))
        stop2("Entry ", i, " of `locusAttributes` has no names")
      if(!all(nms %in% attribNames))
        stop2("Entry ", i, " of `locusAttributes` has illegal entries: ", setdiff(nms, attribNames))
    }
    res = a
  }
  else if (is.list(a) && any(attribNames %in% names(a))) {
    # Format 2: Single list of attributes
    nms = names(a)
    if(!all(nms %in% attribNames))
      stop2("Illegal locus attribute: ", setdiff(nms, attribNames))

    res = list(a)
  }
  else if(is.data.frame(a) || is.matrix(a)) {
    # Format 3: Allelic ladder as data.frame or matrix
    res = freqDb2attribList(a, format = "allelicLadder")
  }
  else if(is.list(a) && all(sapply(a, is.numeric))) {
    # Format 4: Frequency database as a list of frequency vectors
    res = freqDb2attribList(a, format = "list")
  }
  else
    stop2("Unknown format of `locusAttributes")

  res
}

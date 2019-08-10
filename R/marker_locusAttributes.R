#' Get or set locus attributes
#'
#' Retrieve or modify the attributes of attached markers
#'
#' The default setting `markers = NULL` select markers automatically, depending
#' on the `matchNames` argument. If `matchNames = FALSE`, all markers are chosen
#' (and an error issued if this is incompatible with the given
#' `locusAttributes`). If `matchNames = TRUE`, markers will be matched against
#' 'name' attributese included in `locusAttributes` (and an error issued if
#' these are missing).
#'
#' Note that the default value `NA` of `matchNames` is changed to TRUE if all
#' entries of `locusAttributes` have a 'name' component which matches the name a
#' an attached marker.
#'
#' Possible attributes given in `locusAttributes` are as follows (default values
#' in parenthesis):
#'
#' * 'alleles' : a character vector with allele labels
#'
#' * 'afreq' :  a numeric vector with allele frequencies (`rep.int(1/L, L)`,
#' where `L = length(alleles)`)
#'
#' * 'name' : marker name (NA)
#'
#' * 'chrom' : chromosome number (NA)
#'
#' * 'posMb' : physical location in megabases (NA)
#'
#' * 'posCm' : position in centiMorgan (NA)
#'
#' * 'mutmod' : mutation model, or model name (NULL)
#'
#' * 'rate' : mutation model parameter (NULL)
#'
#' @param x A `ped` object, or a list of such.
#' @param markers A character vector (with marker names) or a numeric vector
#'   (with marker indices). If NULL (default), the behaviour depends on
#'   `matchNames`, see Details.
#' @param attribs A subset of the character vector `c("alleles", "afreq", "name"
#'   ,"chrom" ,"posMb","posCm", "mutmod")`.
#' @param locusAttributes A list of lists, with attributes for each marker. If
#'   `locusAttributes` is just a list (not a list of lists), then it is repeated
#'   to match the number of markers.
#' @param matchNames A logical, only relevant if `markers = NULL`. If TRUE, then
#'   the markers to be modified are identified by the 'name' component of each
#'   `locusAttributes` entry. If FALSE, all markers attached to `x` are selected
#'   in order.
#' @param erase A logical. If TRUE, all previous attributes of the selected
#'   markers are erased. If FALSE, attributes not affected by the submitted
#'   `locusAttributes` remain untouched.
#'
#' @return
#'
#' * `getLocusAttributes` : a list of lists
#'
#' * `setLocusAttributes` : a modified version of `x`.
#'
#' @examples
#' x = singleton(1)
#' x = addMarkers(x, marker(x, name = "m1", alleles = 1:2))
#' x = addMarkers(x, marker(x, name = "m2", alleles = letters[1:2], chrom = "X"))
#'
#' # Change frequencies at both loci
#' y = setLocusAttributes(x, markers = 1:2, loc = list(afreq = c(.1, .9)))
#' getMarkers(y, 1)
#'
#' # Set the same mutation model at both loci
#' z = setLocusAttributes(x, markers = 1:2, loc = list(mutmod = "proportional", rate = .1))
#' mutmod(z, 1)
#'
#' # By default, the markers to be modified are identified by name
#' locs = list(list(name = "m1", alleles = 1:10),
#'             list(name = "m2", alleles = letters[1:10]))
#' w = setLocusAttributes(x, loc = locs)
#' getMarkers(w, 1:2)
#'
#' # If `erase = T` attributes not explicitly given are erased
#' w2 = setLocusAttributes(x, loc = locs, erase = TRUE)
#' chrom(w2, 2) # not "X" anymore
#'
#' # The getter and setter are inverses
#' newx = setLocusAttributes(x, loc = getLocusAttributes(x))
#' stopifnot(identical(x, newx))
#'
#' @name locusAttributes
NULL


#' @rdname locusAttributes
#' @export
getLocusAttributes = function(x, markers = NULL,
                              attribs = c("alleles", "afreq", "name", "chrom",
                                          "posMb", "posCm", "mutmod")) {

  if(is.pedList(x))
    x = x[[1]]

  if(!is.ped(x))
    stop2("Input must be a `ped` object or a list of such")

  if(is.null(markers))
    markers = seq_len(nMarkers(x))

  attribs = match.arg(attribs, several.ok = T)

  mlist = getMarkers(x, markers)
  lapply(mlist, function(m) {
    a = attributes(m)[attribs]
    a = a[!is.na(names(a))]
    a
  })
}

#' @rdname locusAttributes
#' @importFrom utils modifyList
#' @export
setLocusAttributes = function(x, markers = NULL, locusAttributes,
                              matchNames = NA, erase = F) {

  # If pedlist input, recurse over components
  if(is.pedList(x)) {
    y = lapply(x, setLocusAttributes, markers = markers,
               locusAttributes = locusAttributes, matchNames = matchNames, erase = erase)
    return(y)
  }

  ### Single `ped` input

  if(!is.ped(x))
    stop2("Input must be a `ped` object or a list of such")

  N = nMarkers(x)
  if(N == 0)
    stop2("This function can only modify already attached markers.\nUse `setMarkers() to attach new markers.")

  # Recycle `locusAttributes` if given as a single list
  recyclingNeeded = is.list(locusAttributes) && !is.list(locusAttributes[[1]])
  if(recyclingNeeded) {
    if(is.null(markers))
      stop2("When `locusAttributes` is a single list, then `markers` cannot be NULL")
    locusAttributes = rep(list(locusAttributes), length(markers))
  }

  # Automatic marker selection
  if(is.null(markers)) {

    if(is.na(matchNames) || isTRUE(matchNames)) {

      # Check if attributes include marker names
      hasNames = all(vapply(locusAttributes, function(a) 'name' %in% names(a), FUN.VALUE = F))
      if(hasNames)
        nms = vapply(locusAttributes, function(a) a[['name']], FUN.VALUE = "")
      else
        nms = names(locusAttributes)

      if(dup <- anyDuplicated(nms))
        stop2("Duplicated marker name in attribute list: ", nms[dup])

      # If matchNames = NA, change to TRUE if all new names match existing ones
      if(is.na(matchNames)) {
        matchNames = !is.null(nms) && all(nms %in% name(x, 1:N))
      }
    }

    # By now, matchNames is either T or F
    if(matchNames) markers = nms
    else markers = 1:N
  }

  if(anyDuplicated(markers))
    stop2("Duplicated markers: ", markers[duplicated(markers)])

  # Index of selected markers
  midx = whichMarkers(x, markers)
  M = length(midx)
  L = length(locusAttributes)

  if(L != M)
    stop2("List of locus attributes does not match the number of markers")

  als = getAlleles(x, markers = midx)
  oldAttrs = getLocusAttributes(x, markers = midx)

  for(i in seq_along(midx)) {
    # Alleles
    ali = als[, c(2*i - 1, 2*i), drop = FALSE]

    # Attributes
    newattri = locusAttributes[[i]]

    if(!erase) {
      updatedattri = modifyList(oldAttrs[[i]], newattri)

      # If new alleles are given without frequencies, the old freqs must be erased anyway
      if("alleles" %in% names(newattri) && !"afreq" %in% names(newattri))
        updatedattri$afreq = NULL

      newattri = updatedattri
    }

    # Create the new marker object (this catches errors!)
    arglist = c(list(x = x, allelematrix = ali), newattri)
    newM = do.call(marker, arglist)

    # Insert in place
    x$markerdata[[midx[i]]] = newM
  }

  # Return modified ped oject
  x
}

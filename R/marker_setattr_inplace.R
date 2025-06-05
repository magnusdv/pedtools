#' Set marker attributes
#'
#' These S3 methods perform in-place modifications of marker attributes. They
#' work on single marker objects and markers attached to ped objects (or lists
#' of such). Although these functions will continue to exist, we recommend the
#' newer alternatives [setGenotype()], [setAfreq()], ... in most cases.
#'
#' @param x Either a `marker` object, a `ped` object or a list of `ped` objects.
#' @param marker,markers The index or name of a marker (or a vector indicating
#'   several markers) attached to `ped`. Used if `x` is a `ped` object.
#' @param id The ID label of a single pedigree member.
#' @param ... Further arguments, not used.
#' @param value Replacement value(s).
#'
#' @return These functions perform in-place modification of `x`.
#'
#' @seealso Alternative setters (not in-place): [marker_setattr].
#'   Marker attribute getters: [marker_getattr].
#'
#' @examples
#' x = nuclearPed(1)
#' x = addMarker(x, alleles = 1:2)
#'
#' # Set genotypes
#' genotype(x, marker = 1, id = 1) = "1/2"
#'
#' # Set marker name
#' name(x, 1) = "M"
#'
#' # Change allele freqs
#' afreq(x, "M") = c(`1` = 0.1, `2` = 0.9)
#'
#' # Set position
#' chrom(x, "M") = 1
#' posMb(x, "M") = 123.45
#'
#' # Check result
#' m = marker(x, `1` = "1/2", name = "M", afreq = c(`1` = 0.1, `2` = 0.9),
#'            chrom = 1, posMb = 123.45)
#' stopifnot(identical(x$MARKERS[[1]], m))
#'
#' @name marker_inplace
NULL



# Genotype ----------------------------------------------------------------

#' @rdname marker_inplace
#' @export
`genotype<-` = function(x, ..., value) {
  UseMethod("genotype<-")
}

#' @rdname marker_inplace
#' @export
`genotype<-.marker` = function(x, id, ..., value) {
  pedlabels = attr(x, 'pedmembers')
  id_int = match(id, pedlabels)

  if (anyNA(id_int))
    stop2("Unknown ID label: ", setdiff(id, pedlabels))

  if(length(value) == 1)
    value = strsplit(as.character(value), "/", fixed = TRUE)[[1]]

  if (!length(value) %in% 1:2)
    stop2("Number of alleles must be 1 or 2: ", value)

  if(length(value) == 1)
    value = rep(value, 2)

  if(length(id_int) > 1)
    value = rep(value, each = length(id_int))

  a = alleles(x)
  g_num = match(value, a, nomatch = 0)

  miss = value[g_num == 0]
  unknown = setdiff(miss, c("0", "", "-", NA))
  if(length(unknown) > 0)
    stop2("Unknown allele for this marker: ", unknown)

  x[id_int, ] = g_num
  x
}

#' @rdname marker_inplace
#' @export
`genotype<-.ped` = function(x, marker, id, ..., value) {
  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Genotype replacement can only be done for a single marker")

  idx = whichMarkers(x, markers = marker)
  genotype(x$MARKERS[[idx]], id) = value
  x
}


# Mutation model ----------------------------------------------------------

#' @rdname marker_inplace
#' @export
`mutmod<-` = function(x, ..., value) {
  UseMethod("mutmod<-")
}

#' @rdname marker_inplace
#' @export
`mutmod<-.marker` = function(x, ..., value) {
  if (!requireNamespace("pedmut", quietly = TRUE))
    stop2("Package `pedmut` must be installed in order to include mutation models")

  if(is.null(value))
    mdl = NULL
  else if(inherits(value, "mutationModel"))
    mdl = value
  else if(is.list(value)) {
    value$alleles = alleles(x)
    value$afreq = afreq(x)
    mdl = do.call(pedmut::mutationModel, value)
  }
  else
    stop2("Mutation model replacement must be either:\n",
          "  * NULL\n",
          "  * a `mutationModel` object\n",
          "  * a list of arguments to be passed onto `pedmut::mutationModel()`")
  attr(x, 'mutmod') = mdl

  x
}

#' @rdname marker_inplace
#' @export
`mutmod<-.ped` = function(x, marker = NULL, ..., value) {
  marker = marker %||% seq_markers(x)

  idx = whichMarkers(x, markers = marker)
  for(i in idx)
    mutmod(x$MARKERS[[i]]) = value

  x
}

#' @rdname marker_inplace
#' @export
`mutmod<-.list` = function(x, marker = NULL, ..., value) {
  for(i in seq_along(x))
    mutmod(x[[i]], marker) = value
  x
}



# Allele frequencies ------------------------------------------------------

#' @rdname marker_inplace
#' @export
`afreq<-` = function(x, ..., value) {
  UseMethod("afreq<-")
}

#' @rdname marker_inplace
#' @export
`afreq<-.marker` = function(x, ..., value) {
  als = alleles(x)
  freqnames = names(value)

  if(is.null(freqnames))
    stop2("Frequency vector must be named (with allele labels)")

  if(anyDuplicated.default(freqnames))
    stop2("Duplicated alleles in frequency vector: ", value[duplicated(value)])

  if(length(freqnames) != length(als) || anyNA(als_order <- match(freqnames, als))) {
    if(length(unknown <- setdiff(freqnames, als)) > 0)
      stop2("Unknown allele: ", unknown)

    if(length(miss <- setdiff(als, freqnames)) > 0)
      stop2("Alleles missing from frequency vector: ", miss)
  }

  if(round(sum(value), 3) != 1)
    stop2("Frequencies must sum to 1")

  attr(x, 'afreq') = as.numeric(value[als_order])
  x
}

#' @rdname marker_inplace
#' @export
`afreq<-.ped` = function(x, marker, ..., value) {
  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Frequency replacement can only be done for a single marker")

  idx = whichMarkers(x, markers = marker)
  afreq(x$MARKERS[[idx]]) = value
  x
}


#' @rdname marker_inplace
#' @export
`afreq<-.list` = function(x, marker, ..., value) {
  for(i in seq_along(x))
    afreq(x[[i]], marker) = value
  x
}

######################
# simple accessors:
# * name()
# * chrom()
# * posMb()
######################


### name setter
#' @rdname marker_inplace
#' @export
`name<-` = function(x, ..., value) {
  UseMethod("name<-")
}

#' @rdname marker_inplace
#' @export
`name<-.marker` = function(x, ..., value) {
  name = as.character(value)

  if(length(name) != 1)
    stop2("Length of `name` must be 1: ", name)
  if (isTRUE(suppressWarnings(name == as.integer(name))))
    stop2("Attribute `name` cannot consist entirely of digits: ", name)

  attr(x, 'name') = name
  x
}

#' @rdname marker_inplace
#' @export
`name<-.ped` = function(x, markers = NULL, ..., value) {
  markers = markers %||% seq_markers(x)
  if(length(markers) == 0)
    stop2("Argument `markers` cannot be empty")
  if(length(value) != length(markers))
    stop2("Length of replacement vector must equal the number of markers")
  if(!is.character(value))
    stop2("Replacement must be a character vector")

  idx = whichMarkers(x, markers = markers)

  # Set names
  x$MARKERS[idx] = lapply(seq_along(idx), function(i) {
    m = x$MARKERS[[idx[i]]]
    name(m) = value[i]
    m
  })

  # Check for duplicates
  checkDupNames(x)

  x
}

#' @rdname marker_inplace
#' @export
`name<-.list` = function(x, markers = NULL, ..., value) {
  lapply(x, function(cmp) `name<-`(cmp, markers = markers, value = value))
}



# Chromosome --------------------------------------------------------------


#' @rdname marker_inplace
#' @export
`chrom<-` = function(x, ..., value) {
  UseMethod("chrom<-")
}

#' @rdname marker_inplace
#' @export
`chrom<-.marker` = function(x, ..., value) {
  chrom = as.character(value)

  if(length(chrom) != 1)
    stop2("Length of `chrom` must be 1: ", chrom)

  attr(x, 'chrom') = chrom
  x
}

#' @rdname marker_inplace
#' @export
`chrom<-.ped` = function(x, markers = NULL, ..., value) {
  markers = markers %||% seq_markers(x)
  if(length(markers) == 0)
    stop2("Argument `markers` cannot be empty")
  if(length(value) > length(markers))
    stop2("Replacement vector larger than the number of markers: ", value)

  if(length(value) < length(markers))
    value = rep_len(value, length(markers))

  idx = whichMarkers(x, markers = markers)

  x$MARKERS[idx] = lapply(seq_along(idx), function(i) {
    m = x$MARKERS[[idx[i]]]
    chrom(m) = value[i]
    m
  })
  x
}

#' @rdname marker_inplace
#' @export
`chrom<-.list` = function(x, markers = NULL, ..., value) {
  lapply(x, function(cmp) `chrom<-`(cmp, markers = markers, value = value))
}



# Position ----------------------------------------------------------------


#' @rdname marker_inplace
#' @export
`posMb<-` = function(x, ..., value) {
  UseMethod("posMb<-")
}


#' @rdname marker_inplace
#' @export
`posMb<-.marker` = function(x, ..., value) {
  pos = suppressWarnings(as.numeric(value))

  if(length(pos) != 1 || (!is.na(value) && is.na(pos)) || (!is.na(pos) && pos < 0))
    stop2("`posMb` replacement must be a single nonnegative number, or NA: ", value)

  attr(x, 'posMb') = pos
  x
}

#' @rdname marker_inplace
#' @export
`posMb<-.ped` = function(x, markers = NULL, ..., value) {
  markers = markers %||% seq_markers(x)

  if(length(markers) == 0)
    stop2("Argument `markers` cannot be empty")

  nm = length(markers)
  nv = length(value)
  if(nv > nm)
    stop2("Replacement vector is longer than the number of markers")
  else if(nv < nm)
    value = rep(value, length.out = nm)

  idx = whichMarkers(x, markers = markers)

  x$MARKERS[idx] = lapply(seq_along(idx), function(i) {
    m = x$MARKERS[[idx[i]]]
    posMb(m) = value[i]
    m
  })
  x
}

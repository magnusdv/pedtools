#' Get/set marker attributes
#'
#' These functions can be used to manipulate a single attribute of one or
#' several markers. Each getter/setter can be used in two ways: Either directly
#' on a `marker` object, or on a `ped` object which has markers attached to it.
#'
#' @param x A [ped] object or a [marker] object
#' @param marker,markers The index or name of a marker (or a vector indicating
#'   several markers) attached to `ped`. Used if `x` is a `ped` object
#' @param id The ID label of a single pedigree member
#' @param ... Further arguments, not used in most of these functions
#' @param value Replacement value(s)
#'
#' @return The getters return the value of the query. The setters perform
#'   in-place modification of the input.
#'
#' @examples
#' x = nuclearPed(1)
#' x = setMarkers(x, locusAttributes = list(name = "M", alleles = 1:2))
#'
#' # Set genotype
#' genotype(x, marker = "M", id = 1) = 1:2
#' genotype(x, marker = "M", id = 3) = 1
#'
#' # Genotypes are returned as a vector of length 2
#' genotype(x, marker = "M", id = 1)
#'
#' # Change allele freqs
#' afreq(x, "M") = c(`1` = 0.1, `2` = 0.9)
#'
#' # Check the new frequencies
#' afreq(x, "M")
#'
#' @name marker_getset
NULL

#' @rdname marker_getset
#' @export
genotype = function(x, ...) {
  UseMethod("genotype")
}

#' @rdname marker_getset
#' @export
genotype.marker = function(x, id, ...) {
  if(length(id) != 1)
    stop2("Argument `id` must have length 1")

  id_int = match(id, attr(x, 'pedmembers'))
  if(is.na(id_int))
    stop2("Unknown ID label: ", id)

  g_num = x[id_int, ]

  g = rep(NA_character_, 2)
  g[g_num > 0] = alleles(x)[g_num]
  g
}

#' @rdname marker_getset
#' @export
genotype.ped = function(x, markers = NULL, id, ...) {
  mlist = getMarkers(x, markers = markers)
  if(length(mlist) == 0)
    stop2("No markers selected")
  if(length(mlist) > 1)
    stop2("More than one marker selected")

  m = mlist[[1]]
  genotype(m, id)
}

#' @rdname marker_getset
#' @export
`genotype<-` = function(x, ..., value) {
  UseMethod("genotype<-")
}

#' @rdname marker_getset
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

#' @rdname marker_getset
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

### mutmod ###
#' @rdname marker_getset
#' @export
mutmod = function(x, ...) {
  UseMethod("mutmod")
}

#' @rdname marker_getset
#' @export
mutmod.marker = function(x, ...) {
  attr(x, 'mutmod')
}

#' @rdname marker_getset
#' @export
mutmod.ped = function(x, marker, ...) {
  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Mutation model can only be accessed for one marker at a time")

  mlist = getMarkers(x, markers = marker)
  m = mlist[[1]]
  mutmod(m)
}

#' @rdname marker_getset
#' @export
mutmod.list = function(x, marker, ...) {
  mutmod(x[[1]], marker = marker)
}

#' @rdname marker_getset
#' @export
`mutmod<-` = function(x, ..., value) {
  UseMethod("mutmod<-")
}

#' @rdname marker_getset
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

#' @rdname marker_getset
#' @export
`mutmod<-.ped` = function(x, marker = NULL, ..., value) {
  marker = marker %||% seq_markers(x)

  idx = whichMarkers(x, markers = marker)
  for(i in idx)
    mutmod(x$MARKERS[[i]]) = value

  x
}

#' @rdname marker_getset
#' @export
`mutmod<-.list` = function(x, marker = NULL, ..., value) {
  marker = marker %||% seq_markers(x)

  for(i in seq_along(x))
    mutmod(x[[i]], marker) = value
  x
}

### alleles ###
#' @rdname marker_getset
#' @export
alleles = function(x, ...) {
  UseMethod("alleles")
}

#' @rdname marker_getset
#' @export
alleles.marker = function(x, ...) {
  attr(x, 'alleles')
}

#' @rdname marker_getset
#' @export
alleles.ped = function(x, marker, ...) {
  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Allele extraction can only be done for a single marker")

  mlist = getMarkers(x, markers = marker)
  m = mlist[[1]]
  alleles(m)
}


### afreq ###
#' @rdname marker_getset
#' @export
afreq = function(x, ...) {
  UseMethod("afreq")
}

#' @rdname marker_getset
#' @export
afreq.marker = function(x, ...) {
  afr = attr(x, "afreq")
  names(afr) = alleles(x)
  afr
}

#' @rdname marker_getset
#' @export
afreq.ped = function(x, marker, ...) {
  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Frequency extraction can only be done for a single marker")
  mlist = getMarkers(x, markers = marker)

  m = mlist[[1]]
  afreq(m)
}

#' @rdname marker_getset
#' @export
afreq.list = function(x, marker, ...) {
  comp_wise = lapply(x, afreq, marker = marker)
  if(!listIdentical(comp_wise))
    stop2("The output of `afreq()` differs between pedigree components")
  comp_wise[[1]]
}

#' @rdname marker_getset
#' @export
`afreq<-` = function(x, ..., value) {
  UseMethod("afreq<-")
}

#' @rdname marker_getset
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

  attr(x, 'afreq') = value[als_order]
  x
}

#' @rdname marker_getset
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


#' @rdname marker_getset
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


### getters
#' @rdname marker_getset
#' @export
name = function(x, ...) {
  UseMethod("name")
}

#' @rdname marker_getset
#' @export
name.marker = function(x, ...) {
  attr(x, 'name')
}

#' @rdname marker_getset
#' @export
name.ped = function(x, markers, ...) {
  mlist = getMarkers(x, markers = markers)
  vapply(mlist, name.marker, character(1))
}

#' @rdname marker_getset
#' @export
name.list = function(x, markers, ...) {
  comp_wise = lapply(x, name, markers = markers)
  if(!listIdentical(comp_wise))
    stop2("The output of `name()` differs between pedigree components")
  comp_wise[[1]]
}

### name setter
#' @rdname marker_getset
#' @export
`name<-` = function(x, ..., value) {
  UseMethod("name<-")
}

#' @rdname marker_getset
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

#' @rdname marker_getset
#' @export
`name<-.ped` = function(x, markers, ..., value) {
  if(missing(markers) || length(markers) == 0)
    stop2("Argument `markers` cannot be empty")
  if(length(value) != length(markers))
    stop2("Length of replacement vector must equal the number of markers")
  if(!is.character(value))
    stop2("Replacement must be a character vector")
  if(anyDuplicated.default(value))
    stop2("Replacement values must be unique")

  idx = whichMarkers(x, markers = markers)

  x$MARKERS[idx] = lapply(seq_along(idx), function(i) {
    m = x$MARKERS[[idx[i]]]
    name(m) = value[i]
    m
  })
  x
}

#' @rdname marker_getset
#' @export
`name<-.list` = function(x, markers, ..., value) {
  lapply(x, function(cmp) `name<-`(cmp, markers = markers, value = value))
}


#' @rdname marker_getset
#' @export
chrom = function(x, ...) {
  UseMethod("chrom")
}

#' @rdname marker_getset
#' @export
chrom.marker = function(x, ...) {
  attr(x, 'chrom')
}

#' @rdname marker_getset
#' @export
chrom.ped = function(x, markers, ...) {
  mlist = getMarkers(x, markers = markers)
  vapply(mlist, chrom.marker, character(1))
}

#' @rdname marker_getset
#' @export
chrom.list = function(x, markers, ...) {
  comp_wise = lapply(x, chrom, markers = markers)
  if(!listIdentical(comp_wise))
    stop2("The output of `chrom()` differs between pedigree components")
  comp_wise[[1]]
}

### chrom setter
#' @rdname marker_getset
#' @export
`chrom<-` = function(x, ..., value) {
  UseMethod("chrom<-")
}

#' @rdname marker_getset
#' @export
`chrom<-.marker` = function(x, ..., value) {
  chrom = as.character(value)

  if(length(chrom) != 1)
    stop2("Length of `chrom` must be 1: ", chrom)

  attr(x, 'chrom') = chrom
  x
}

#' @rdname marker_getset
#' @export
`chrom<-.ped` = function(x, markers, ..., value) {
  if(missing(markers) || length(markers) == 0)
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

#' @rdname marker_getset
#' @export
`chrom<-.list` = function(x, markers, ..., value) {
  lapply(x, function(cmp) `chrom<-`(cmp, markers = markers, value = value))
}

#' @rdname marker_getset
#' @export
posMb = function(x, ...) {
  UseMethod("posMb")
}

#' @rdname marker_getset
#' @export
posMb.marker = function(x, ...) {
  as.numeric(attr(x, 'posMb'))
}

#' @rdname marker_getset
#' @export
posMb.ped = function(x, markers, ...) {
  mlist = getMarkers(x, markers = markers)
  vapply(mlist, posMb, numeric(1))
}


### posMb setter
#' @rdname marker_getset
#' @export
`posMb<-` = function(x, ..., value) {
  UseMethod("posMb<-")
}

#' @rdname marker_getset
#' @export
`posMb<-.marker` = function(x, ..., value) {
  pos = suppressWarnings(as.numeric(value))

  if((!is.na(value) && is.na(pos)) || length(pos) != 1)
    stop2("`posMb` replacement must be a single number: ", value)
  if(pos < 0)
    stop2("`posMb` replacement must be nonnegative: ", value)

  attr(x, 'posMb') = pos
  x
}

#' @rdname marker_getset
#' @export
`posMb<-.ped` = function(x, markers, ..., value) {
  if(missing(markers) || length(markers) == 0)
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

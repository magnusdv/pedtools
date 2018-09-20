### genotype ###

#' @export
genotype = function(x, ...) {
  UseMethod("genotype")
}

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

#' @export
genotype.ped = function(x, markers=NULL, id, ...) {
  mlist = getMarkers(x, markers=markers)
  if(length(mlist) == 0)
    stop2("No markers selected") #TODO
  if(length(mlist) > 1)
    stop2("More than one marker selected")

  m = mlist[[1]]
  genotype(m, id)
}

#' @export
`genotype<-` = function(x, ..., value) {
  UseMethod("genotype<-")
}

#' @export
`genotype<-.marker` = function(x, id, ..., value) {
  pedlabels = attr(x, 'pedmembers')
  id_int = match(id, pedlabels)

  if (anyNA(id_int))
    stop2("Unknown ID label: ", setdiff(id, pedlabels))

  if (!length(value) %in% 1:2)
    stop2("Length of genotype vector must be 1 or 2")

  if(length(value) == 1)
    value = rep(value, 2)

  if(length(id_int) > 1)
    value = rep(value, each=length(id_int))

  a = alleles(x)
  g_num = match(value, a, nomatch=0)

  miss = value[g_num == 0]
  unknown = setdiff(miss, c("0", "", "-", NA))
  if(length(unknown) > 0)
    stop2("Unknown allele for this marker: ", unknown)

  x[id_int, ] = g_num
  x
}

#' @export
`genotype<-.ped` = function(x, marker, id, ..., value) {
  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Genotype replacement can only be done for a single marker")

  idx = whichMarkers(x, markers=marker)
  genotype(x$markerdata[[idx]], id) = value
  x
}

### mutmod ###
#' @export
mutmod = function(x, ...) {
  UseMethod("mutmod")
}

#' @export
mutmod.marker = function(x, ...) {
  attr(x, 'mutmod')
}

#' @export
mutmod.ped = function(x, marker, ...) {
  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Mutation model can only be accessed for one marker at a time")

  mlist = getMarkers(x, markers=marker)
  m = mlist[[1]]
  mutmod(m)
}

#' @export
`mutmod<-` = function(x, ..., value) {
  UseMethod("mutmod<-")
}

#' @export
`mutmod<-.marker` = function(x, ..., value) {
  if (!requireNamespace("pedmut", quietly = TRUE))
    stop2("Package `pedmut` must be installed in order to include mutation models")

  als = alleles(x)
  attr(x, 'mutmod') = pedmut::mutationModel(model = value, alleles = als)
  x
}

#' @export
`mutmod<-.ped` = function(x, marker, ..., value) {
  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Mutation model replacement can only be done for a single marker")

  idx = whichMarkers(x, markers=marker)
  mutmod(x$markerdata[[idx]]) = value
  x
}

### alleles ###

#' @export
alleles = function(x, ...) {
  UseMethod("alleles")
}

#' @export
alleles.marker = function(x, ...) {
  attr(x, 'alleles')
}

#' @export
alleles.ped = function(x, marker, ...) {
  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Allele extraction can only be done for a single marker")

  mlist = getMarkers(x, markers=marker)
  m = mlist[[1]]
  alleles(m)
}


### afreq ###

#' @export
afreq = function(x, ...) {
  UseMethod("afreq")
}

#' @export
afreq.marker = function(x, ...) {
  afr = attr(x, "afreq")
  names(afr) = alleles(x)
  afr
}

#' @export
afreq.ped = function(x, marker, ...) {
  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Frequency extraction can only be done for a single marker")
  mlist = getMarkers(x, markers=marker)

  m = mlist[[1]]
  afreq(m)
}

#' @export
`afreq<-` = function(x, ..., value) {
  UseMethod("afreq<-")
}

#' @export
`afreq<-.marker` = function(x, ..., value) {
  als = alleles(x)
  freqnames = names(value)

  if(is.null(freqnames))
    stop2("Frequency vector must be named (with allele labels)")

  if(anyDuplicated(freqnames))
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

#' @export
`afreq<-.ped` = function(x, marker, ..., value) {
  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Frequency replacement can only be done for a single marker")

  idx = whichMarkers(x, markers=marker)
  afreq(x$markerdata[[idx]]) = value
  x
}

######################
# simple accessors:
# * name()
# * chrom()
# * posMb()
# * posCm()
######################


### getters

#' @export
name = function(x, ...) {
  UseMethod("name")
}

#' @export
name.marker = function(x, ...) {
  attr(x, 'name')
}

#' @export
name.ped = function(x, markers, ...) {
  mlist = getMarkers(x, markers=markers)
  vapply(mlist, name, character(1))
}

#' @export
chrom = function(x, ...) {
  UseMethod("chrom")
}

#' @export
chrom.marker = function(x, ...) {
  attr(x, 'chrom')
}

#' @export
chrom.ped = function(x, markers, ...) {
  mlist = getMarkers(x, markers=markers)
  vapply(mlist, chrom, character(1))
}

#' @export
posMb = function(x, ...) {
  UseMethod("posMb")
}

#' @export
posMb.marker = function(x, ...) {
  as.numeric(attr(x, 'posMb'))
}

#' @export
posMb.ped = function(x, markers, ...) {
  mlist = getMarkers(x, markers=markers)
  vapply(mlist, posMb, numeric(1))
}

#' @export
posCm = function(x, ...) {
  UseMethod("posCm")
}

#' @export
posCm.marker = function(x, ...) {
  as.numeric(attr(x, 'posCm'))
}

#' @export
posCm.ped = function(x, markers, ...) {
  mlist = getMarkers(x, markers=markers)
  vapply(mlist, posCm, numeric(1))
}


### name setter
#' @export
`name<-` = function(x, ..., value) {
  UseMethod("name<-")
}

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

#' @export
`name<-.ped` = function(x, markers, ..., value) {
  if(missing(markers) || length(markers) == 0)
    stop2("Argument `markers` cannot be empty")
  if(length(value) != length(markers))
    stop2("Length of replacement vector must equal the number of markers")
  if(!is.character(value))
    stop2("Replacement must be a character vector")
  if(anyDuplicated(value))
    stop2("Replacement values must be unique")

  idx = whichMarkers(x, markers=markers)

  x$markerdata[idx] = lapply(seq_along(idx), function(i) {
    m = x$markerdata[[idx[i]]]
    name(m) = value[i]
    m
  })
  x
}

### chrom setter
#' @export
`chrom<-` = function(x, ..., value) {
  UseMethod("chrom<-")
}

#' @export
`chrom<-.marker` = function(x, ..., value) {
  chrom = as.character(value)

  if(length(chrom) != 1)
    stop2("Length of `chrom` must be 1: ", chrom)

  attr(x, 'chrom') = chrom
  x
}

#' @export
`chrom<-.ped` = function(x, markers, ..., value) {
  if(missing(markers) || length(markers) == 0)
    stop2("Argument `markers` cannot be empty")
  if(length(value) != length(markers))
    stop2("Length of replacement vector must equal the number of markers")
  if(anyDuplicated(value))
    stop2("Replacement values must be unique")

  idx = whichMarkers(x, markers=markers)

  x$markerdata[idx] = lapply(seq_along(idx), function(i) {
    m = x$markerdata[[idx[i]]]
    chrom(m) = value[i]
    m
  })
  x
}

### posCm setter

#' @export
`posCm<-` = function(x, ..., value) {
  UseMethod("posCm<-")
}

#' @export
`posCm<-.marker` = function(x, ..., value) {
  pos = suppressWarnings(as.numeric(value))

  if((!is.na(value) && is.na(pos)) || length(pos) != 1)
    stop2("`posCm` replacement must be a single number: ", value)
  if(pos < 0)
    stop2("`posCm` replacement must be nonnegative: ", value)

  attr(x, 'posCm') = pos
  x
}

#' @export
`posCm<-.ped` = function(x, markers, ..., value) {
  if(missing(markers) || length(markers) == 0)
    stop2("Argument `markers` cannot be empty")
  if(length(value) != length(markers))
    stop2("Length of replacement vector must equal the number of markers")
  if(anyDuplicated(value))
    stop2("Replacement values must be unique")

  idx = whichMarkers(x, markers=markers)

  x$markerdata[idx] = lapply(seq_along(idx), function(i) {
    m = x$markerdata[[idx[i]]]
    posCm(m) = value[i]
    m
  })
  x
}

### posMb setter

#' @export
`posMb<-` = function(x, ..., value) {
  UseMethod("posMb<-")
}

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

#' @export
`posMb<-.ped` = function(x, markers, ..., value) {
  if(missing(markers) || length(markers) == 0)
    stop2("Argument `markers` cannot be empty")
  if(length(value) != length(markers))
    stop2("Length of replacement vector must equal the number of markers")
  if(anyDuplicated(value))
    stop2("Replacement values must be unique")

  idx = whichMarkers(x, markers=markers)

  x$markerdata[idx] = lapply(seq_along(idx), function(i) {
    m = x$markerdata[[idx[i]]]
    posMb(m) = value[i]
    m
  })
  x
}

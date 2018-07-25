### genotype ###

#' @export
genotype = function(x, ...) {
  UseMethod("genotype")
}

#' @export
genotype.marker = function(x, id, ...) {
  if(length(id) != 1)
    stop("Argument `id` must have length 1", call.=F)

  id_int = match(id, attr(x, 'pedmembers'))
  if(is.na(id_int))
    stop("Unknown ID label: ", id, call.=F)

  g_num = x[id_int, ]

  g = rep(NA_character_, 2)
  g[g_num > 0] = alleles(x)[g_num]
  g
}

#' @export
genotype.ped = function(x, markers=NULL, id, ...) {
  mlist = getMarkers(x, markers=markers)
  if(length(mlist) == 0)
    stop("No markers selected") #TODO
  if(length(mlist) > 1)
    stop("More than one marker selected")

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
    stop("Unknown ID label: ", paste(setdiff(id, pedlabels), collapse=", "),
         call.=F)

  if (!length(value) %in% 1:2)
    stop("Length of genotype vector must be 1 or 2", call.=F)

  if(length(value) == 1)
    value = rep(value, 2)

  if(length(id_int) > 1)
    value = rep(value, each=length(id_int))

  a = alleles(x)
  g_num = match(value, a, nomatch=0)

  miss = value[g_num == 0]
  unknown = setdiff(miss, c("0", "", "-", NA))
  if(length(unknown) > 0)
    stop("Unknown allele for this marker: ", paste(unknown, collapse=", "),
         call.=F)

  x[id_int, ] = g_num
  x
}

#' @export
`genotype<-.ped` = function(x, marker, id, ..., value) {
  if(missing(marker) || length(marker) == 0)
    stop("Argument `marker` cannot be empty", call.=F)
  if(length(marker) > 1)
    stop("Genotype replacement can only be done for a single marker", call.=F)

  idx = whichMarkers(x, markers=marker)
  genotype(x$markerdata[[idx]], id) = value
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
    stop("Argument `marker` cannot be empty", call.=F)
  if(length(marker) > 1)
    stop("Allele extraction can only be done for a single marker", call.=F)

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
    stop("Argument `marker` cannot be empty", call.=F)
  if(length(marker) > 1)
    stop("Frequency extraction can only be done for a single marker", call.=F)
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
    stop("Frequency vector must be named (with allele labels)", call.=F)

  if(anyDuplicated(freqnames))
    stop("Duplicated alleles in frequency vector: ", value[duplicated(value)], call.=F)

  if(length(freqnames) != length(als) || anyNA(als_order <- match(freqnames, als))) {
    if(length(unknown <- setdiff(freqnames, als)) > 0)
      stop("Unknown allele: ", paste(unknown, collapse=", "), call.=F)

    if(length(miss <- setdiff(als, freqnames)) > 0)
      stop("Alleles missing from frequency vector: ", paste(miss, collapse=", "), call.=F)
  }

  if(round(sum(value), 3) != 1)
    stop("Frequencies must sum to 1", call.=F)

  attr(x, 'afreq') = value[als_order]
  x
}

#' @export
`afreq<-.ped` = function(x, marker, ..., value) {
  if(missing(marker) || length(marker) == 0)
    stop("Argument `marker` cannot be empty", call.=F)
  if(length(marker) > 1)
    stop("Frequency replacement can only be done for a single marker", call.=F)

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
  as.character(attr(x, 'name'))
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
  as.character(attr(x, 'chrom'))
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


### setters
#' @export
`name<-` = function(x, ..., value) {
  UseMethod("name<-")
}

#' @export
`name<-.marker` = function(x, ..., value) {
  if(!is.character(value) || length(value) != 1)
    stop("Replacement value must be a string", call.=F)
  attr(x, 'name') = value
  x
}

#' @export
`name<-.ped` = function(x, markers, ..., value) {
  if(missing(markers) || length(markers) == 0)
    stop("Argument `marker` cannot be empty", call.=F)
  if(length(value) != length(markers))
    stop("Length of replacement vector must equal the number of markers", call.=F)
  if(!is.character(value))
    stop("Replacement must be a character vector", call.=F)
  if(anyDuplicated(value))
    stop("Replacement values must be unique", call.=F)

  idx = whichMarkers(x, markers=markers)

  x$markerdata[idx] = lapply(seq_along(idx), function(i) {
    m = x$markerdata[[idx[i]]]
    name(m) = value[i]
    m
  })
  x
}


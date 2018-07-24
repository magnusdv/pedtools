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
    stop("ID label not associated with this marker: ", id, call.=F)

  g_num = x[id_int, ]

  g = rep(NA_character_, 2)
  g[g_num > 0] = alleles(x)[g_num]
  g
}

#' @export
genotype.ped = function(x, markeridx=NULL, markername=NULL, id, ...) {
  mlist = getMarkers(x, markeridx=markeridx, markername=markername)
  if(length(mlist) == 0)
    stop("No markers selected")
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
  id_int = match(id, attr(x, 'pedmembers'))

  if(anyNA(id_int))
    stop("ID label not associated with this marker: ",
         paste(setdiff(id, x$LABELS), collapse=", "), call.=F)

  if(length(value)==1)
    value = rep(value, 2)

  if(length(id_int)==1 && length(value) != 2)
    stop("When replacing the genotype of a single individual, ",
         "the replacement vector must have length either 1 or 2", call.=F)

  if(length(id_int) > 1)
    if(length(value) == 2)
      value = rep(value, each=length(id_int))
    else if(length(value) != 2*length(id_int))
      stop("Wrong length of replacement vector. See ?genotype", call.=F)

    a = alleles(x)
    g_num = match(value, a, nomatch=0)

    miss = value[g_num == 0]
    if(!all(miss %in% c("0", "", "-", NA)))
      stop("Alleles associated with this marker are: ", paste(a, collapse=", "), call.=F)

    x[id_int, ] = g_num
    x
}

#' @export
`genotype<-.ped` = function(x, markeridx=NULL, markername=NULL, id, ..., value) {
  if(is.null(markeridx))
    markeridx = whichMarkers(x, markername=markername)
  if(length(markeridx) == 0)
    stop("No markers selected", call.=F)
  if(length(markeridx) > 1)
    stop("Genotype replacement can only be done for a single marker", call.=F)

  genotype(x$markerdata[[markeridx]], id) = value
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
alleles.ped = function(x, markeridx=NULL, markername=NULL, ...) {
  mlist = getMarkers(x, markeridx=markeridx, markername=markername)
  if(length(mlist) == 0)
    stop("No markers selected")
  if(length(mlist) > 1)
    stop("More than one marker selected")

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
  attr(x, "afreq")
}

#' @export
afreq.ped = function(x, markeridx=NULL, markername=NULL, ...) {
  mlist = getMarkers(x, markeridx=markeridx, markername=markername)
  if(length(mlist) == 0)
    stop("No markers selected")
  if(length(mlist) > 1)
    stop("More than one marker selected")

  m = mlist[[1]]
  afreq(m)
}

#' @export
`afreq<-` = function(x, ..., value) {
  UseMethod("afreq<-")
}

#' @export
`afreq<-.marker` = function(x, ..., value) {
  old = attr(x, 'afreq')
  if(length(value) != length(old))
    stop("Replacement vector must have same length as the original: ",
         length(old), call.=F)
  if(round(sum(value), 3) != 1)
    stop("Frequency vector must sum to 1", call.=F)
  attr(x, 'afreq') = value
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
name.ped = function(x, markeridx=NULL, ...) {
  mlist = getMarkers(x, markeridx=markeridx)
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
chrom.ped = function(x, markeridx=NULL, markername=NULL, ...) {
  mlist = getMarkers(x, markeridx=markeridx, markername=markername)
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
posMb.ped = function(x, markeridx=NULL, markername=NULL, ...) {
  mlist = getMarkers(x, markeridx=markeridx, markername=markername)
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
posCm.ped = function(x, markeridx=NULL, markername=NULL, ...) {
  mlist = getMarkers(x, markeridx=markeridx, markername=markername)
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
`name<-.ped` = function(x, markeridx=NULL, markername=NULL, ..., value) {
  if(!is.character(value))
    stop("Replacement must be a character vector", call.=F)

  if(is.null(markeridx))
    markeridx = whichMarkers(x, markername=markername)
  if(length(markeridx) == 0)
    stop("No markers selected", call.=F)
  if(length(value) != length(markeridx))
    stop("Length of replacement vector must equal the number of markers", call.=F)
  if(anyDuplicated(value))
    stop("Replacement values must be unique", call.=F)

  x$markerdata[markeridx] = lapply(seq_along(markeridx), function(i) {
    m = x$markerdata[[markeridx[i]]]
    name(m) = value[i]
    m
  })
  x
}


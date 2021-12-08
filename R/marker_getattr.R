#' Get marker attributes
#'
#' S3 methods retrieving marker attributes. They work on single marker objects
#' and markers attached to ped objects (or lists of such).
#'
#' @param x Either a `marker` object, a `ped` object or a list of `ped` objects.
#' @param marker,markers The index or name of a marker (or a vector indicating
#'   several markers) attached to `x`.
#' @param id The ID label of a single pedigree member.
#' @param ... Further arguments, not used.
#'
#' @return The associated marker attributes.
#'
#' @seealso Setting marker attributes: [marker_setattr] and [marker_inplace].
#'
#' @examples
#' x = nuclearPed(1)
#' x = addMarker(x) # add empty marker
#'
#' # Inspect default attributes
#' alleles(x, marker = 1)
#' afreq(x, marker = 1)
#' name(x, marker = 1)  # NA
#' chrom(x, marker = 1) # NA
#'
#' @name marker_getattr
NULL


# Get genotype ------------------------------------------------------------


#' @rdname marker_getattr
#' @export
genotype = function(x, ...) {
  UseMethod("genotype")
}

#' @rdname marker_getattr
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

#' @rdname marker_getattr
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


# NOT SURE ABOUT THIS YET
getGenotypeMarker = function(m, id, sep = "/") {
  pedlabels = attr(m, 'pedmembers')
  id_int = match(id, pedlabels)

  if (anyNA(id_int))
    stop2("Unknown ID label: ", setdiff(id, pedlabels))

  als = alleles(m)

  a1int =  m[id_int, 1]
  a1int[a1int == 0L] = NA
  a1 = als[a1int]

  a2int =  m[id_int, 2]
  a2int[a2int == 0L] = NA
  a2 = als[a2int]

  if(is.null(sep)) {
    if(length(id) != 1)
      return(c(a1, a2))
    else
      stop2("When `sep` is NULL, `id` must have length 1")
  }

  res = paste(a1, a2, sep = sep)

  res
}

# NOT SURE ABOUT THIS YET
getGenotype = function(x, marker = NULL, id, sep = "/") {
  if(is.null(marker))
    marker = seq_len(nMarkers(x))

  if(length(marker) == 0)
    stop2("No marker identified")

  nma = length(marker)
  nid = length(id)
  if(nma != 1 && nid != 1)
    stop2("Either `marker` or `id` must have length 1")

  if(is.null(sep) && any(c(nma, nid) > 1))
    stop2("When `sep` is NULL, both `marker` and `id` must have length 1")

  if(nid == 1 && is.pedList(x)) {
    comp = getComponent(x, id, checkUnique = TRUE, errorIfUnknown = TRUE)
    x = x[[comp]]
  }

  ### End of prep

  if(nid == 1 && nma == 1) {
    m = getMarkers(x, markers = marker)[[1]]
    return(getGenotypeMarker(m, id, sep = sep))
  }

  if(nid == 1) {
    mlist = getMarkers(x, markers = marker)
    return(vapply(mlist, function(m) getGenotypeMarker(m, id, sep = sep), "1"))
  }

  if(nma == 1) {
    if(is.pedList(x)) {
      comp = getComponent(x, id, checkUnique = TRUE, errorIfUnknown = TRUE)
      res = vapply(seq_along(id), function(k)
        getGenotype(x[[comp[k]]], marker = marker, id = id[k], sep = sep), "1")
      return(res)
    }

    m = getMarkers(x, markers = marker)[[1]]
    return(getGenotypeMarker(m, id = id, sep = sep))
  }
}



# Get mutation model ------------------------------------------------------

#' @rdname marker_getattr
#' @export
mutmod = function(x, ...) {
  UseMethod("mutmod")
}

#' @rdname marker_getattr
#' @export
mutmod.marker = function(x, ...) {
  attr(x, 'mutmod')
}

#' @rdname marker_getattr
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

#' @rdname marker_getattr
#' @export
mutmod.list = function(x, marker, ...) {
  mutmod(x[[1]], marker = marker)
}



# Get alleles -------------------------------------------------------------

#' @rdname marker_getattr
#' @export
alleles = function(x, ...) {
  UseMethod("alleles")
}

#' @rdname marker_getattr
#' @export
alleles.marker = function(x, ...) {
  attr(x, 'alleles')
}

#' @rdname marker_getattr
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

#' @rdname marker_getattr
#' @export
alleles.list = function(x, marker, ...) {
  comp_wise = lapply(x, alleles, marker = marker)
  if(!listIdentical(comp_wise))
    stop2("The output of `alleles()` differs between pedigree components")
  comp_wise[[1]]
}


# Get allele frequencies --------------------------------------------------

#' @rdname marker_getattr
#' @export
afreq = function(x, ...) {
  UseMethod("afreq")
}

#' @rdname marker_getattr
#' @export
afreq.marker = function(x, ...) {
  afr = attr(x, "afreq")
  names(afr) = alleles(x)
  afr
}

#' @rdname marker_getattr
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

#' @rdname marker_getattr
#' @export
afreq.list = function(x, marker, ...) {
  comp_wise = lapply(x, afreq, marker = marker)
  if(!listIdentical(comp_wise))
    stop2("The output of `afreq()` differs between pedigree components")
  comp_wise[[1]]
}


# Get marker name ---------------------------------------------------------

#' @rdname marker_getattr
#' @export
name = function(x, ...) {
  UseMethod("name")
}

#' @rdname marker_getattr
#' @export
name.marker = function(x, ...) {
  attr(x, 'name')
}

#' @rdname marker_getattr
#' @export
name.ped = function(x, markers = NULL, ...) {
  markers = markers %||% seq_markers(x)

  mlist = getMarkers(x, markers = markers)
  vapply(mlist, name.marker, character(1))
}

#' @rdname marker_getattr
#' @export
name.list = function(x, markers = NULL, ...) {
  comp_wise = lapply(x, name, markers = markers)
  if(!listIdentical(comp_wise))
    stop2("The output of `name()` differs between pedigree components")
  comp_wise[[1]]
}


# Get chromosome ----------------------------------------------------------

#' @rdname marker_getattr
#' @export
chrom = function(x, ...) {
  UseMethod("chrom")
}

#' @rdname marker_getattr
#' @export
chrom.marker = function(x, ...) {
  attr(x, 'chrom')
}

#' @rdname marker_getattr
#' @export
chrom.ped = function(x, markers = NULL, ...) {
  markers = markers %||% seq_markers(x)

  mlist = getMarkers(x, markers = markers)
  vapply(mlist, chrom.marker, character(1))
}

#' @rdname marker_getattr
#' @export
chrom.list = function(x, markers = NULL, ...) {
  comp_wise = lapply(x, chrom, markers = markers)
  if(!listIdentical(comp_wise))
    stop2("The output of `chrom()` differs between pedigree components")
  comp_wise[[1]]
}


# Position ----------------------------------------------------------------

#' @rdname marker_getattr
#' @export
posMb = function(x, ...) {
  UseMethod("posMb")
}

#' @rdname marker_getattr
#' @export
posMb.marker = function(x, ...) {
  as.numeric(attr(x, 'posMb'))
}

#' @rdname marker_getattr
#' @export
posMb.ped = function(x, markers = NULL, ...) {
  markers = markers %||% seq_markers(x)

  mlist = getMarkers(x, markers = markers)
  vapply(mlist, posMb, numeric(1))
}

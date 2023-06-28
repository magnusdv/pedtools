#' Set marker attributes
#'
#' These functions set or modify various attributes of markers attached to a
#' pedigree. They are sometimes more convenient (and pipe-friendly) than the
#' in-place modifiers described in [marker_inplace].
#'
#' @inheritParams marker
#' @param x A `ped` object or a list of `ped` objects.
#' @param id,ids A vector naming one or several pedigree members, or a function
#'   (e.g., [founders()]).
#' @param marker A vector of indices or names of one or several markers attached
#'   to `x`.
#' @param name A character of the same length as `marker`, containing marker
#'   names.
#' @param chrom A character of the same length as `marker`, containing
#'   chromosome labels.
#' @param posMb A numeric of the same length as `marker`, containing the
#'   physical marker positions in megabases (or NA).
#' @param strict A logical. If TRUE (default) the new frequencies cannot remove
#'   or add any alleles.
#'
#' @return A copy of `x` with modified attributes.
#'
#' @examples
#' x = nuclearPed(1) |>
#'   addMarker(alleles = 1:2) |>
#'   setMarkername(marker = 1, name = "M") |>
#'   setGenotype(marker = "M", ids = 1, geno = "1/2") |>
#'   setAfreq(marker = "M", afreq = c(`1` = 0.1, `2` = 0.9)) |>
#'   setChrom(marker = "M", chrom = 1) |>
#'   setPosition(marker = "M", posMb = 123.45)
#'
#' # Of course, all of this could have been done on creation:
#' y = addMarker(nuclearPed(),
#'               `1` = "1/2",
#'               afreq = c(`1` = 0.1, `2` = 0.9),
#'               name = "M",
#'               chrom = 1,
#'               posMb = 123.45)
#' stopifnot(identical(x, y))
#'
#' @name marker_setattr
NULL

# Set genotype ------------------------------------------------------------

setGenotypeMarker = function(m, ids, geno) {
  pedlabels = attr(m, 'pedmembers')
  id_int = match(ids, pedlabels)

  if (anyNA(id_int))
    stop2("Unknown ID label: ", setdiff(ids, pedlabels))

  nid = length(ids)
  als = alleles(m)

  # Special case: geno gives two separate alleles
  # Exception: X marker and two males
  if(length(geno) == 2 && all(geno %in% als) && !(isXmarker.marker(m) && nid == 2 && all(attr(m, "sex")[id_int] == 1)))
    geno = paste(geno, collapse = "/")

  ng = length(geno)

  if(nid > 1 && ng == 1)
    geno = rep(geno, nid)
  else if(nid != ng)
    stop2("When setting the genotype of multiple individuals, `geno` must have length either 1 or `length(ids)`")

  genoList = strsplit(as.character(geno), "/", fixed = TRUE)

  lg = lengths(genoList)

  # Error if no alleles
  empt = lg == 0
  if(any(empt))
    stop2("No alleles given for individual ", ids[which(empt)])

  # Error if too many alleles
  bad = lg > 2
  if(any(bad)) {
    idx = which(bad)[1]
    stop2(sprintf("Too many alleles given for individual '%s': ", ids[idx]), genoList[[idx]])
  }

  # Error if only 1 allele given (except X males)
  short = lg == 1
  if(isXmarker.marker(m))
    short[attr(m, "sex")[id_int] == 1] = FALSE
  if(any(short))
    stop2("Only one allele given for individual ", ids[which(short)])

  # Remaining short (X males!): Duplicate
  genoList[lg == 1] = lapply(genoList[lg == 1], rep_len, 2)

  gAls = unlist(genoList, use.names = FALSE)
  gInt = match(gAls, als, nomatch = 0)

  # Check for unknown alleles
  miss = gAls[gInt == 0]
  unknown = setdiff(miss, c("0", "", "-", NA))
  if(length(unknown) > 0) {
    nm = attr(m, "name")
    stop2(sprintf("Unknown allele for %s: ", ifelse(is.na(nm), "this marker", nm)) , unknown)
  }

  # Insert alleles
  even = seq_len(nid) * 2
  m[id_int, 1] = gInt[even - 1]
  m[id_int, 2] = gInt[even]

  m
}


#' @rdname marker_setattr
#' @export
setGenotype = function(x, marker = NULL, ids = NULL, geno = NULL, id = NULL) {
  if(!is.null(id))
    ids = id

  if(is.null(marker))
    marker = seq_len(nMarkers(x))

  nma = length(marker)

  if(is.function(ids))
    ids = ids(x)

  nid = length(ids)
  if(nma != 1 && nid != 1)
    stop2("Either `marker` or `ids` must have length 1")

  if(nma == 1 && nid == 1) {

    if(is.pedList(x)) {
      comp = getComponent(x, ids, checkUnique = TRUE, errorIfUnknown = TRUE)
      x[[comp]] = setGenotype(x[[comp]], marker = marker, ids = ids, geno = geno)
      return(x)
    }

    midx = whichMarkers(x, markers = marker)
    x$MARKERS[[midx]] = setGenotypeMarker(x$MARKERS[[midx]], ids = ids, geno = geno)
    return(x)
  }

  if(nma == 1 ) { # length(ids) > 1
    if(length(geno) == 1)
      geno = rep(geno, nid)
    else if(length(geno) != nid)
      stop("When setting the genotype of multiple individuals, `geno` must have length either 1 or `length(ids)`")

    if(is.pedList(x)) {
      comp = getComponent(x, ids, checkUnique = TRUE, errorIfUnknown = TRUE)
      x[[comp]] = lapply(x[[comp]], function(y) {
        idx = match(ids, y$ID, nomatch = 0)
        setGenotype(y, marker = marker, ids = ids[idx], geno = geno[idx])
      })

      return(x)
    }

    midx = whichMarkers(x, markers = marker)
    x$MARKERS[[midx]] = setGenotypeMarker(x$MARKERS[[midx]], ids = ids, geno = geno)
    return(x)
  }

  if(nid == 1) { # length(marker) > 1

    if(length(geno) == 1)
      geno = rep(geno, nma)
    else if(length(geno) != nma)
      stop("When setting the genotype for multiple markers, `geno` must have length either 1 or `length(marker)`")

    if(is.pedList(x)) {
      comp = getComponent(x, ids, checkUnique = TRUE, errorIfUnknown = TRUE)
      x[[comp]] = setGenotype(x[[comp]], marker = marker, ids = ids, geno = geno)
      return(x)
    }

    midx = whichMarkers(x, markers = marker)
    x$MARKERS[midx] = lapply(seq_along(midx), function(k) {
      m = x$MARKERS[[midx[k]]]
      setGenotypeMarker(m, ids = ids, geno = geno[k])
    })
    return(x)
  }

}


# Set allele frequencies --------------------------------------------------

#' @rdname marker_setattr
#' @export
setAfreq = function(x, marker, afreq, strict = TRUE) {

  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Frequency replacement can only be done for a single marker")

  if(is.pedList(x)) { # leads to lots of redundant tests, but no big deal
    y = lapply(x, setAfreq, marker = marker, afreq = afreq)
    return(y)
  }

  idx = whichMarkers(x, markers = marker)
  m = x$MARKERS[[idx]]
  als = alleles(m)

  freqnames = names(afreq) %||% stop2("The frequency vector must be named with allele labels")

  if(anyDuplicated.default(freqnames))
    stop2("Duplicated alleles in frequency vector: ", afreq[duplicated(freqnames)])

  if(round(sum(afreq), 3) != 1)
    stop2("Frequencies must sum to 1")

  if(strict) {
    if(length(miss <- setdiff(als, freqnames)))
      stop2("Alleles missing from frequency vector: ", miss)

    alsOrder = match(freqnames, als)
    if(anyNA(alsOrder))
      stop2("Unknown allele: ", freqnames[is.na(alsOrder)])
    afreq = afreq[alsOrder]
    attr(m, "afreq") = as.numeric(afreq)
  }
  else {
    if(allowsMutations(m))
      warning("Mutation models may be invalidated by frequency change")

    # Sort (numerically if appropriate)
    if (!anyNA(suppressWarnings(as.numeric(freqnames))))
      alsOrder = order(as.numeric(freqnames))
    else
      alsOrder = order(freqnames)
    afreq = afreq[alsOrder]
    freqnames = names(afreq)

    # Recompute index matrix
    obs = als[m[m > 0]] # observed alleles
    m[m > 0] = match(obs, freqnames)
    if(anyNA(m))
      stop2("Invalid allele: ", setdiff(obs, freqnames))

    attr(m, "afreq") = unname(afreq)
    attr(m, "alleles") = freqnames
  }

  x$MARKERS[[idx]] = m
  x
}



# Set marker name ---------------------------------------------------------

#' @rdname marker_setattr
#' @export
setMarkername = function(x, marker = NULL, name) {

  if(is.pedList(x)) { # leads to redundant tests, but no big deal
    y = lapply(x, setMarkername, marker = marker, name = name)
    return(y)
  }

  marker = marker %||% seq_markers(x)
  if(length(marker) == 0)
    stop2("Argument `marker` cannot be empty")

  name = as.character(name) # especially important to allow NA input

  if(length(name) != length(marker))
    stop2("Length of `name` must equal the number of markers")
  if(!is.character(name))
    stop2("Argument `name` must be a character vector")
  if (any(dig <- suppressWarnings(name == as.integer(name)), na.rm = TRUE))
    stop2("Marker name cannot consist entirely of digits: ", name[dig])

  idx = whichMarkers(x, markers = marker)
  mlist = x$MARKERS[idx]
  x$MARKERS[idx] = lapply(seq_along(idx), function(i) `attr<-`(mlist[[i]], "name", name[i]))

  # Check for duplicates
  checkDupNames(x)

  x
}


# Set chromosome ----------------------------------------------------------


#' @rdname marker_setattr
#' @export
setChrom = function(x, marker = NULL, chrom) {

  if(is.pedList(x)) { # leads to redundant tests, but no big deal
    y = lapply(x, setChrom, marker = marker, chrom = chrom)
    return(y)
  }

  marker = marker %||% seq_markers(x)
  if(length(marker) == 0)
    stop2("Argument `marker` cannot be empty")

  chrom = as.character(chrom)
  if(length(chrom) != length(marker))
    if(length(chrom) == 1)
      chrom = rep(chrom, length(marker))
  else
    stop2("Length of `chrom` must equal the number of markers")

  idx = whichMarkers(x, markers = marker)
  mlist = x$MARKERS[idx]
  x$MARKERS[idx] = lapply(seq_along(idx), function(i) `attr<-`(mlist[[i]], "chrom", chrom[i]))

  x
}


# Set position ----------------------------------------------------------------


#' @rdname marker_setattr
#' @export
setPosition = function(x, marker = NULL, posMb) {

  if(is.pedList(x)) { # leads to redundant tests, but no big deal
    y = lapply(x, setPosition, marker = marker, posMb = posMb)
    return(y)
  }

  marker = marker %||% seq_markers(x)
  if(length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(posMb) != length(marker))
    stop2("Length of `posMb` must equal the number of markers")

  pos = suppressWarnings(as.numeric(posMb))
  bad = (!is.na(posMb) & is.na(pos)) | (!is.na(pos) & pos < 0)
  if(any(bad))
    stop2("`posMb` must contain nonnegative numbers or NA: ", posMb[bad])

  idx = whichMarkers(x, markers = marker)
  mlist = x$MARKERS[idx]
  x$MARKERS[idx] = lapply(seq_along(idx), function(i) `attr<-`(mlist[[i]], "posMb", pos[i]))

  x
}


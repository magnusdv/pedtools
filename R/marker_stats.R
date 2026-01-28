#' Find homozygous genotypes
#'
#' Identifies homozygous genotypes in the marker data. A genotype is homozygous
#' if both alleles are non-missing and equal.
#'
#' @param x A `ped` object or a list of such. An error is raised if `x` has no
#'   marker data.
#' @param ids A vector of individual ID labels. Defaults to all typed members of
#'   `x`.
#' @param count A logical. If TRUE, return counts per individual.
#'
#' @returns By default, a logical matrix of dimension `N x L`, where N is
#'   `length(ids)` and L is the number of markers. If `count = TRUE`, a numeric
#'   vector of length N giving the number of homozygous genotypes for each
#'   individual.
#'
#' @examples
#' x = nuclearPed(father = "fa", mother = "mo", children = "ch") |>
#'   addMarker(name = "m1", geno = c(NA, "1/1", "1/2")) |>
#'   addMarker(name = "m2", geno = c(NA, "2/2", "2/2"))
#'
#' isHomozygous(x)
#' isHomozygous(x, ids = "mo")
#' isHomozygous(x, count = TRUE)
#'
#' @export
isHomozygous = function(x, ids = typedMembers(x), count = FALSE) {
  if(!is.ped(x))
    return(.isHomozygous0(x, ids = ids, count = count))

  nMark = nMarkers(x, compwise = FALSE)
  if(nMark == 0) stop2("Input pedigree has no marker data")

  amat = unlist(x$MARKERS)
  dim(amat) = c(length(x$ID), 2 * nMark)

  idsInt = internalID(x, ids)

  odd = seq(1, dim(amat)[2], by = 2)
  even = odd + 1

  a1 = amat[idsInt, odd, drop = FALSE]
  a2 = amat[idsInt, even, drop = FALSE]

  homoz = a1 == a2
  homoz[a1 == 0] = NA  # set missing to NA

  # Fix dim names
  nms = name(x)
  nms[is.na(nms)] = which(is.na(nms))
  colnames(homoz) = nms
  rownames(homoz) = ids

  if(count)
    return(rowSums(homoz, na.rm = TRUE))

  homoz
}

# Old version, slower, but takes ped lists
.isHomozygous0 = function(x, ids = typedMembers(x), count = FALSE) {
  if(!hasMarkers(x)) stop2("Input pedigree has no marker data")

  amat = getAlleles(x, ids)
  ids = rownames(amat)

  odd = seq(1, dim(amat)[2], by = 2)
  even = odd + 1

  homoz = amat[, odd, drop = FALSE] == amat[, even, drop = FALSE]

  # Fix column names
  nms = name(x)
  nms[is.na(nms)] = which(is.na(nms))
  colnames(homoz) = nms

  if(count)
    return(rowSums(homoz, na.rm = TRUE))

  homoz
}


#' Find markers for which two individuals have the same genotype
#'
#' Identifies markers for which two individuals have the same (non-missing)
#' genotype. The comparison is done after sorting the genotypes internally.
#'
#' @param x A `ped` object or a list of such. An error is raised if `x` has no
#'   marker data.
#' @param ids A vector of two individual ID labels.
#' @param count A logical. If TRUE, return the number of markers with shared
#'   genotype.
#'
#' @returns A logical vector with one entry per marker (`NA` if either genotype
#'   is missing). If `count = TRUE`, the number of `TRUE` entries.
#'
#' @seealso [isHomozygous()], [sortGenotypes()]
#'
#' @examples
#' x = nuclearPed() |>
#'   addMarker(name = "m1", geno = c(NA, "1/1", "1/2")) |>
#'   addMarker(name = "m2", geno = c(NA, "1/2", "2/1"))
#'
#' sameGenotype(x, 2:3)
#' sameGenotype(x, 2:3, count = TRUE)
#'
#' @export
sameGenotype = function(x, ids = typedMembers(x), count = FALSE) {
  if(!is.ped(x))
    return(.sameGenotype0(x, ids = ids, count = count))

  nMark = nMarkers(x, compwise = FALSE)
  if(nMark == 0) stop2("Input pedigree has no marker data")
  if(length(ids) != 2) stop2("Argument `ids` must be a vector of length 2")

  amat = unlist(x$MARKERS)
  dim(amat) = c(length(x$ID), 2 * nMark)

  idsInt = internalID(x, ids)
  id1 = idsInt[1]
  id2 = idsInt[2]

  odd = seq(1, dim(amat)[2], by = 2)
  even = odd + 1

  a1 = amat[id1, odd]
  a2 = amat[id1, even]
  b1 = amat[id2, odd]
  b2 = amat[id2, even]

  sameG = (a1 == b1 & a2 == b2) | (a1 == b2 & a2 == b1)
  sameG[a1+a2 == 0] = NA  # set missing to NA

  # Fix dim names
  nms = name(x)
  if(anyNA(nms))
    nms[is.na(nms)] = which(is.na(nms))
  names(sameG) = nms

  if(count)
    return(sum(sameG, na.rm = TRUE))

  sameG
}


# Old version, slower, but takes ped lists
.sameGenotype0 = function(x, ids = typedMembers(x), count = FALSE) {
  if(!hasMarkers(x)) stop2("Input pedigree has no marker data")
  if(length(ids) != 2) stop2("Argument `ids` must be a vector of length 2")

  x = sortGenotypes(x)

  gmat = getGenotypes(x, ids)
  same = gmat[1, ] == gmat[2, ]

  if(count)
    return(sum(same, na.rm = TRUE))

  same
}


#' Expected homozygosity and heterozygosity
#'
#' Computes the expected homozygosity and heterozygosity for a marker from its
#' allele frequencies. The homozygosity is \eqn{\sum_i p_i^2} and the
#' heterozygosity is one minus this.
#'
#' @param p A numeric vector of allele frequencies, or a `marker` object (from
#'   which the frequencies are extracted with [afreq()]).
#'
#' @returns A numeric value giving the expected homozygosity or heterozygosity.
#'
#' @examples
#' p = c(0.2, 0.5, 0.3)
#' expectedHomozygosity(p)
#' expectedHeterozygosity(p)
#'
#' @export
expectedHomozygosity = function(p) {
  if(is.marker(p))
    p = afreq(p)
  else if(!is.numeric(p) || any(p < 0) || round(sum(p), 3) != 1)
    stop2("`p` is neither a `marker` object nor a valid frequency vector: ", p)
  sum(p^2)
}

#' @rdname expectedHomozygosity
#' @export
expectedHeterozygosity = function(p) {
  if(is.marker(p))
    p = afreq(p)
  else if(!is.numeric(p) || any(p < 0) || round(sum(p), 3) != 1)
    stop2("`p` is neither a `marker` object nor a valid frequency vector: ", p)
  1 - sum(p^2)
}

#' Count homozygous markers
#'
#' Counts the number of homozygous genotypes for each individual.
#'
#' @param x A pedigree (ped object) with marker data.
#' @param ids A vector of individual ID labels. If NULL (default) all
#'   individuals are included.
#'
#' @returns A named integer vector with the counts of homozygous markers for
#'   each individual in `ids`.
#'
#' @examples
#' x = nuclearPed() |> addMarker(`2` = "1/1", `3` = "1/2")
#' countHomozygous(x)
#'
#' @export
countHomozygous = function(x, ids = NULL) {
  if(!hasMarkers(x)) stop2("Input pedigree has no marker data")

  amat = getAlleles(x, ids)
  ids = rownames(amat)

  dm = dim(amat)
  odd = seq(1, dm[2], by = 2)
  even = odd + 1

  homoz = amat[, odd] == amat[, even]
  res = .rowSums(homoz, dm[1], dm[2]/2, na.rm = TRUE)
  names(res) = ids
  res
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

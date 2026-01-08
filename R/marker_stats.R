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

#' Allele frequency database
#'
#' Extract marker allele frequencies as a data.frame, sorted according to the
#' 'allelic ladder'
#'
#' @param x A `ped` object
#' @param markers A character vector (with marker names) or a numeric vector
#'   (with marker indices)
#'
#' @return A data.frame with allele labels as now names, and marker labels as
#'   column names
#'
#' @examples
#' loc1 = list(name = "m1", alleles = 3:4, afreq = c(.1, .9))
#' loc2 = list(name = "m2", alleles = c("1", "10.2", "3"))
#' x = setMarkers(singleton(1), locus = list(loc1, loc2))
#' getFrequencyDatabase(x)
#'
#' @export
getFrequencyDatabase = function(x, markers) {
  if(missing(markers))
    markers = seq_len(nMarkers(x))

  mlist = getMarkers(x, markers)
  mnames = unlist(lapply(mlist, name))
  afr = lapply(mlist, afreq)
  allAlleles = unique.default(names(unlist(afr)))

  nums = !anyNA(suppressWarnings(is.numeric(allAlleles)))
  if(nums)
    allAlleles = allAlleles[order(as.numeric(allAlleles))]
  else
    allAlleles = allAlleles[order(allAlleles)]

  res = matrix(NA, nrow = length(allAlleles), ncol = length(mlist),
               dimnames = list(allAlleles, mnames))

  for(i in seq_along(mlist)) {
    a = afr[[i]]
    res[names(a), i] = a
  }

  res
}


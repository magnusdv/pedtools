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
getFrequencyDatabase = function(x, markers = seq_len(nMarkers(x))) {

  # If list of pedigrees, use the first component
  if(is.pedList(x))
    return(getFrequencyDatabase(x[[1]], markers))

  attrs = getLocusAttributes(x, markers, attribs = c("name", "alleles", "afreq"))
  mnames = lapply(attrs, '[[', "name")
  mnames = unlist(mnames) # what about NA's?
  als = lapply(attrs, '[[', "alleles")
  afr = lapply(attrs, '[[', "afreq")

  allAlleles = unique.default(unlist(als))

  # Sort numerically if possible
  nums = !anyNA(suppressWarnings(is.numeric(allAlleles)))
  if(nums)
    allAlleles = allAlleles[order(as.numeric(allAlleles))]
  else
    allAlleles = allAlleles[order(allAlleles)]

  # Create output matrix
  res = matrix(NA, nrow = length(allAlleles), ncol = length(attrs),
               dimnames = list(allAlleles, mnames))

  # Fill inn correct frequencies for each marker
  for(i in seq_along(attrs)) {
    als.i = als[[i]]
    afr.i = afr[[i]]
    res[als.i, i] = afr.i
  }

  res
}

setFrequencyDatabase = function(x, database) {

  database = as.data.frame(database)
  als = rownames(database)
  mnames = colnames(database)

  loci = lapply(database, function(m) {
    idx = !is.na(m)
    list(alleles = als[idx], afreq = m[idx])
  })

  loci
}

## Internal methods
getLocusAttributes = function(x, markers = seq_len(nMarkers(x)),
                   attribs = c("alleles", "afreq", "name" ,"chrom" ,"posMb","posCm")) {

  attribs = match.arg(attribs, several.ok = T)
  mlist = getMarkers(x, markers)
  lapply(mlist, function(m) attributes(m)[attribs])
}


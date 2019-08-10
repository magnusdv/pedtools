#' Allele frequency database
#'
#' Extract or set marker allele frequencies using the 'allelic ladder' format
#'
#' These functions are essentially special cases of [getLocusAttributes()] and
#' [setLocusAttributes()].
#'
#' If `database` is a data.frame, it must have the allele labels as row names
#' and marker names as column names.
#'
#' If `database` is a file path, it is read with the command
#' `read.table(database, header = T, row.names = 1, as.is = T, check.names = F, ...)`
#'
#' @param x A `ped` object, or a list of such
#' @param markers A character vector (with marker names) or a numeric vector
#'   (with marker indices)
#' @param database Either a data.frame or a character file path. See details.
#' @param ... Optional arguments passed on to [read.table()]
#' @return
#'
#' * `getFrequencyDatabase` : a data.frame with allele labels as row names, and
#' marker labels as column names
#'
#' * `getFrequencyDatabase`: a modified version of `x`
#'
#' @seealso  [setLocusAttributes()], [setMarkers()], [setAlleles()]
#'
#' @examples
#' loc1 = list(name = "m1", alleles = 3:4, afreq = c(.1, .9))
#' loc2 = list(name = "m2", alleles = c("1", "10.2", "3"))
#' x = setMarkers(singleton(1), locus = list(loc1, loc2))
#' db = getFrequencyDatabase(x)
#' db
#'
#' y = setFrequencyDatabase(x, database = db)
#' stopifnot(identical(x, y))
#'
#' # The database can also be read directly from file
#' tmp = tempfile()
#' write.table(db, tmp, row.names = TRUE, col.names = TRUE)
#'
#' z = setFrequencyDatabase(x, database = tmp)
#' stopifnot(all.equal(x, z))
#'
#' @name freqDatabase
NULL


#' @rdname freqDatabase
#' @export
getFrequencyDatabase = function(x, markers = NULL) {

  # If list of pedigrees, use the first component
  if(is.pedList(x))
    return(getFrequencyDatabase(x[[1]], markers = markers))

  if(!is.ped(x))
    stop2("Input must be a `ped` object or a list of such")

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

#' @rdname freqDatabase
#' @export
setFrequencyDatabase = function(x, database, ...) {

  # Read the table if file path
  if(is.character(database) && length(database) == 1) {
    database = read.table(database, header = T, row.names = 1, as.is = T, check.names = F, ...)
  }
  else {
    database = as.data.frame(database)
  }

  als = rownames(database)
  nms = colnames(database)
  if(dup <- anyDuplicated(nms))
    stop2("Duplicated marker name in database file: ", nms[dup])


  loci = lapply(seq_along(database), function(i) {

    # Column of frequencies
    frqs = database[[i]]

    # Index of rows with non-missing entries
    idx = !is.na(frqs)

    list(name = nms[i], alleles = als[idx], afreq = frqs[idx])
  })

  setLocusAttributes(x, locusAttributes = loci, matchNames = T, erase = F)
}


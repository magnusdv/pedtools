#' Allele frequency database
#'
#' Functions for reading, setting and extracting allele frequency databases, in
#' either "list" format or "allelic ladder" format.
#'
#' A frequency database in "list" format is a list of numeric vectors; each
#' vector named with the allele labels, and the list itself named with the
#' marker names.
#'
#' Text files containing frequencies in "list" format should look as follows,
#' where "marker1" and "marker2" are marker names, and "a1","a2",... are allele
#' labels (which may be characters or numeric, but will always be converted to
#' characters):
#' \preformatted{
#' marker1
#' a1 0.2
#' a2 0.5
#' a3 0.3
#'
#' marker2
#' a1 0.9
#' a2 0.1
#' }
#'
#' A database in "allelic ladder" format is rectangular, i.e., a numeric matrix
#' (or data frame), with allele labels as row names and markers as column
#' names. `NA` entries correspond to unobserved alleles.
#'
#' @param x	A ped object, or a list of such.
#' @param markers	A character vector (with marker names) or a numeric vector
#'   (with marker indices).
#' @param database Either a list or matrix/data frame with allele frequencies,
#'   or a file path (to be passed on to `readFreqDatabase()`).
#' @param format Either "list" or "ladder".
#' @param filename The path to a text file containing allele frequencies either
#'   in "list" or "allelic ladder" format.
#' @param ... Optional arguments passed on to [read.table()].
#' @return
#'
#' * `getFreqDatabase`: either a list (if `format = "list"`) or a data
#' frame (if `format = "ladder"`)
#'
#' * `readFreqDatabase`: a list (also if `format = "ladder"`) of named
#' numeric vectors
#'
#' * `setFreqDatabase`: a modified version of `x`
#'
#' @seealso  [setLocusAttributes()], [setMarkers()], [setAlleles()]
#'
#' @examples
#' loc1 = list(name = "m1", afreq = c(a = .1, b = .9))
#' loc2 = list(name = "m2", afreq = c("1" = .2, "10.2" = .3, "3" = .5))
#' x = setMarkers(singleton(1), locus = list(loc1, loc2))
#' db = getFreqDatabase(x)
#' db
#'
#' y = setFreqDatabase(x, database = db)
#' stopifnot(identical(x, y))
#'
#' # The database can also be read directly from file
#' tmp = tempfile()
#' write("m1\na 0.1\nb 0.9\n\nm2\n1 0.2\n3 0.5\n10.2 0.3", tmp)
#'
#' z = setFreqDatabase(x, database = tmp)
#' stopifnot(all.equal(x, z))
#'
#' @name freqDatabase
NULL


#' @rdname freqDatabase
#' @export
getFreqDatabase = function(x, markers = NULL, format = c("list", "ladder")) {

  # If list of pedigrees, use the first component
  if(is.pedList(x))
    return(getFreqDatabase(x[[1]], markers = markers, format = format))

  if(!is.ped(x))
    stop2("Input must be a `ped` object or a list of such")

  format = match.arg(format)

  attrs = getLocusAttributes(x, markers, attribs = c("name", "alleles", "afreq"))
  mnames = lapply(attrs, '[[', "name")
  mnames = unlist(mnames) # what about NA's?

  als = lapply(attrs, '[[', "alleles")
  afr = lapply(attrs, '[[', "afreq")

  if(format == "list") {
    res = lapply(seq_along(attrs), function(i) setNames(afr[[i]], als[[i]]))
    names(res) = mnames
  }
  if(format == "ladder") {
    allAlleles = unique.default(unlist(als))

    # Sort numerically if possible
    if(!anyNA(suppressWarnings(is.numeric(allAlleles))))
      allAlleles = allAlleles[order(as.numeric(allAlleles))]
    else
      allAlleles = allAlleles[order(allAlleles)]

    # Create output matrix
    res = matrix(NA, nrow = length(allAlleles), ncol = length(attrs),
                 dimnames = list(allAlleles, mnames))

    # Fill inn correct frequencies for each marker
    for(i in seq_along(attrs)) {
      res[als[[i]], i] = afr[[i]]
    }
  }

  res
}

#' @rdname freqDatabase
#' @export
setFreqDatabase = function(x, database, format = c("list", "ladder"), ...) {

  if(!hasMarkers(x))
    stop2("This function can only modify already attached markers.",
          "\nUse `setMarkers() to attach new markers.")

  format = match.arg(format)

  # If file path, read database
  if(is.character(database) && length(database) == 1) {
    database = readFreqDatabase(database, format = format, ...)
    format = "list" # output is always in list format
  }

  loci = freqDb2attribList(database, format = format)

  # Check matching
  xMarkers = name(x, seq_markers(x))
  dbMarkers = sapply(loci, '[[', "name")

  mtch = match(xMarkers, dbMarkers, nomatch = 0L)
  if(!any(mtch > 0))
    stop2("No matching entries in database")
  if(!all(mtch > 0))
    message("Skipping marker not found in database: ", toString(xMarkers[mtch == 0]))

  useLoci = loci[mtch]
  setLocusAttributes(x, locusAttributes = useLoci, matchNames = TRUE, erase = FALSE)
}


# Conversion frequency database to list of locus attributes
freqDb2attribList = function(database, format = c("list", "ladder")) {

  format = match.arg(format)
  if(format == "list") {
    attrList = lapply(names(database), function(nm) {
      v = database[[nm]]
      als = names(v)
      if(is.null(als) || !is.numeric(v)) {
        stop2("When `database` is a list, then each entry must:\n",
              " * be a numeric vector of allele frequencies\n",
              " * be named with allelic labels")
      }
      list(name = nm, alleles = als, afreq = unname(v))
    })
  }
  else if(format == "ladder") {
    als = rownames(database)
    nms = colnames(database)
    if(dup <- anyDuplicated(nms))
      stop2("Duplicated marker name in database file: ", nms[dup])

    attrList = lapply(seq_along(nms), function(i) {
      frqs = database[, i]
      idx = !is.na(frqs)  # Index of non-missing entries
      list(name = nms[i], alleles = als[idx], afreq = frqs[idx])
    })
  }

  attrList
}


#' @rdname freqDatabase
#' @export
readFreqDatabase = function(filename, format = c("list", "ladder"), ...) {

  format = match.arg(format)

  if(format == "list") {
    # Read as table with two columns
    raw = read.table(filename, colClasses = c("character", "numeric"),
                     header = FALSE, fill = TRUE, blank.lines.skip = FALSE, ...)
    raw = raw[raw$V1 != "", 1:2, drop = FALSE]

    # First/last lines numbers for each marker
    newMarker = which(is.na(raw$V2))
    stops = c(newMarker[-1] -  1, nrow(raw))
    nms = raw[newMarker, 1]

    # Convert to list of named frequency vectors
    res = lapply(seq_along(nms), function(i) {
      m = raw[newMarker[i]:stops[i], , drop = FALSE]
      als = m[-1, 1]
      afr = m[-1, 2]
      setNames(afr, als)
    })
  }
  else if(format == "ladder") {

    # Read entire data frame
    database = read.table(filename, header = TRUE, row.names = 1,
                          as.is = TRUE, check.names = FALSE, ...)
    als = rownames(database)
    nms = colnames(database)

    res = lapply(seq_along(nms), function(i) {
      frqs = database[, i] # Column of frequencies
      idx = !is.na(frqs)   # Index of rows with non-missing entries
      setNames(frqs[idx], als[idx])
      })
  }

  # Add marker names
  names(res) = nms
  res
}

#' @rdname freqDatabase
#' @export
writeFreqDatabase = function(x, filename, markers = NULL, format = c("list", "ladder")) {
  if(is.ped(x) || is.pedList(x))
    db = getFreqDatabase(x, markers, format = "list")
  else
    db = x

  format = match.arg(format)
  if(format == "ladder")
    stop2("Ladder format not implemented yet")

  N = length(db)
  for(i in seq_len(N)) {
    write(names(db)[[i]], file = filename, append = i > 1)
    m = db[[i]]
    write.table(cbind(names(m), unname(m)), file = filename, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
    if(i < N)
      write("", file = filename, append = TRUE)
  }
}

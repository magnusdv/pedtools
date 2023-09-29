#' Allele frequency database
#'
#' Functions for reading, setting and extracting allele frequency databases, in
#' either "list" format, "merlin" format or "allelic ladder" format.
#'
#' A frequency database in "list" format is a list of numeric vectors; each
#' vector named with the allele labels, and the list itself named with the
#' marker names.
#'
#' Text files containing frequencies in "list" format should look as follows,
#' where "M1" and "M2" are marker names, and "a1","a2",... are allele labels
#' (which may be characters or numeric, but will always be converted to
#' characters):
#' \preformatted{
#' M1
#' a1 0.2
#' a2 0.5
#' a3 0.3
#'
#' M2
#' a1 0.9
#' a2 0.1
#' }
#'
#' In "merlin" format, used by the software MERLIN (Abecasis et. al, 2002), the
#' same frequency data would be presented as follows:
#' \preformatted{
#' M M1
#' A a1 0.2
#' A a2 0.5
#' A a3 0.3
#' M M2
#' A a1 0.9
#' A a2 0.1
#' }
#'
#' A database in "allelic ladder" format is rectangular, i.e., a numeric matrix
#' (or data frame), with allele labels as row names and markers as column names.
#' `NA` entries correspond to unobserved alleles.
#'
#' @param x	A ped object, or a list of such.
#' @param markers	A character vector (with marker names) or a numeric vector
#'   (with marker indices).
#' @param database Either a list or matrix/data frame with allele frequencies,
#'   or a file path (to be passed on to `readFreqDatabase()`).
#' @param format Either "list", "ladder" or "merlin" (only in
#'   `readFreqDatabase()`).
#' @param fixNames A logical, by default FALSE. If TRUE all marker names are
#'   converted to upper case, and all periods and space characters are replaced
#'   with "_" (underscore).
#' @param scale1 A logical, by default FALSE. If TRUE, all frequency vectors are
#'   scaled to ensure that it sums to 1.
#' @param filename The path to a text file containing allele frequencies either
#'   in "list" or "allelic ladder" format.
#' @param df A data frame of allele frequencies in either "list" or "allelic
#'   ladder" format. This can be supplied instead of `filename`.
#' @param verbose A logical.
#' @param ... Optional arguments passed on to [read.table()], e.g. `sep = "\t"`
#'   if the file is tab separated.
#'
#' @return
#'
#' * `getFreqDatabase`: either a list (if `format = "list"`) or a data
#' frame (if `format = "ladder"`).
#'
#' * `readFreqDatabase`: a list of named numeric vectors.
#'
#' * `setFreqDatabase`: a modified version of `x`.
#'
#' @seealso  [setLocusAttributes()], [setMarkers()], [setAlleles()].
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

    nms = names(database)
    if(is.null(nms))
      stop2("Frequency database does not include marker names")
    if(dup <- anyDuplicated(nms))
      stop2("Duplicated marker name in frequency database: ", nms[dup])

    attrList = lapply(nms, function(nm) {
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
readFreqDatabase = function(filename = NULL, df = NULL, format = c("list", "ladder", "merlin"),
                            fixNames = FALSE, scale1 = FALSE, verbose = TRUE, ...) {

  format = match.arg(format)

  if(verbose) {
    if(is.null(df))
      cat("Reading frequency file:", filename, "\n")
    cat("Database format:", format, "\n")
  }

  db = switch(format,
    list = { # data frame with 2 columns

      if(is.null(df)) {
        df = read.table(filename, colClasses = c("character", "numeric"),
                        header = FALSE, fill = TRUE, blank.lines.skip = FALSE, ...)
      }

      # Skip blank lines
      blank = df[[1]] == "" | is.na(df[[1]])
      df = df[!blank, 1:2, drop = FALSE]

      # First/last lines numbers for each marker
      labs = df[[1]]
      freqs = df[[2]]

      # New marker starts where freq column is empty (1st col = marker name)
      nameRow = which(is.na(freqs))
      nms = labs[nameRow]

      freqStart = nameRow + 1
      freqStop = c(nameRow[-1] -  1, length(freqs))

      # Convert to list of named frequency vectors
      lapply(seq_along(nms), function(i) {
        rws = freqStart[i]:freqStop[i]
        setNames(freqs[rws], labs[rws])
      })
    },
    merlin = { # Similar to above, but 3 columns; first column contains "M", "A"

      if(is.null(df)) {
        df = read.table(filename, colClasses = c("character", "character", "numeric"),
                        header = FALSE, fill = TRUE, ...)
      }

      # New marker starts where first column has "M"
      nameRow = which(df$V1 == "M")
      nms = df$V2[nameRow]

      freqStart = nameRow + 1
      freqStop = c(nameRow[-1] -  1, nrow(df))

      # Convert to list of named frequency vectors
      lapply(seq_along(nms), function(i) {
        rws = freqStart[i]:freqStop[i]
        setNames(df$V3[rws], df$V2[rws])
      })
    },
    ladder = {  # allelic ladder

      if(is.null(df)) {
        df = read.table(filename, header = TRUE, row.names = 1, as.is = TRUE,
                        check.names = FALSE, ...)
      }
      mat = as.matrix(df)
      als = rownames(mat)
      nms = colnames(mat)

      lapply(seq_along(nms), function(i) {
        frqs = mat[, i]    # column of frequencies
        idx = !is.na(frqs) # rows with non-missing entries
        setNames(frqs[idx], als[idx])
      })
    })

  if(fixNames) {
    # Default fixes: (i) convert all to upper case, (ii) replace space and "." with "_"
    if(verbose)
      cat("Fixing names:\n")

    oldnms = nms
    nms = sub(" ", "_", sub(".", "_", toupper(oldnms), fixed = TRUE), fixed = TRUE)
    if(verbose) {
      changed = which(nms != oldnms)
      for(i in changed)
        cat(sprintf(" %s -> %s\n", oldnms[i], nms[i]))
      if(!length(changed))
        cat(" No changes done")
    }
  }

  # Add marker names
  names(db) = nms

  # Check frequencies
  db = lapply(nms, function(m) {
    afr = db[[m]]

    num = suppressWarnings(as.numeric(afr))
    if(anyNA(num))
      stop2(sprintf("Non-numeric frequency detected for `%s`: %s", m, afr[is.na(num)]))
    afr[] = num

    # Fix bad microvariant labels, e.g. "9.30000000007"
    numals = suppressWarnings(as.numeric(names(afr)))
    if(!anyNA(numals))
      names(afr) = as.character(numals)

    # Scale to sum 1
    if(scale1) {
      sm = sum(afr)
      if(sm != 1) {
        if(verbose) cat(sprintf("Scaling `%s` to 1. Original sum: %g\n", m, sm))
        afr = afr/sm
      }
    }

    afr
  })

  # Add marker names (again...)
  names(db) = nms

  db
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

library(stringr)

#' Read a pedigree from a set of files
#'
#' Assumptions:
#' * All markers present in .dat, .freq and .map files are also present in the .ped file and in the same order.
#' * The .dat file only contains a list of markers and is thus only used for validation.
#'
#' @param prefix a filepath without extension from which .ped, .map, .freq and
#'   .dat files will be read.
#' @param what a subset of the vector c("ped", "map", "dat", "freq") indicating
#'   which files to read. By default all files that are present will be read.
#' @param verbose wheteher to show descriptive messages while reading the
#'   pedigree. Defaults to TRUE.
#' @return a `ped` object
#'
#' @author Elias Hernandis <eliashernandis@gmail.com>
#'
#' @seealso [writePed()]
#' @export

readPed = function(prefix, what = NULL, verbose = TRUE) {
  ped = NULL
  pedFileName = paste(prefix, '.ped', sep = '')
  if ("ped" %in% what || is.null(what)) {
    if (verbose) print(sprintf("Reading %s...", pedFileName))

    df = read.table(pedFileName)
    colnames(df) = c('famid', 'id', 'fid', 'mid', 'sex')
    newdf = df[,1:5]

    # read genotype data if available
    if (ncol(df) >= 6) {
      for (i in 6:ncol(df)) {
        splt = gsub('-', '0', stringr::str_split_fixed(df[,i], '/', 2))
        newdf[6 + 2 * (i - 6)] = splt[,1]
        newdf[6 + 2 * (i - 6) + 1] = splt[,2]
      }
    }
    ped = as.ped(newdf)
  } else {
    stop2("Cannot load pedigree without .ped file")
  }

  mapFileName = paste(prefix, '.map', sep = '')
  if ('map' %in% what || (is.null(what) && file.exists(mapFileName) && file.size(mapFileName) > 0)) {
    if (!file.exists(mapFileName)) {
      stop2("Could not find requested .map file")
    }

    if (verbose) print(sprintf("Reading %s...", mapFileName))
    df = read.table(mapFileName, stringsAsFactors = FALSE)

    for (i in 1:nrow(df)) {
      attr(ped$markerdata[[i]], 'chrom')  = df[i,1]
      attr(ped$markerdata[[i]], 'name')   = df[i,2]
      # TODO: import marker position
    }
  }

  freqFileName = paste(prefix, '.freq', sep = '')
  if ('freq' %in% what || (is.null(what) && file.exists(freqFileName) && file.size(freqFileName) > 0)) {
    if (!file.exists(freqFileName)) {
      stop2("Could not find requested .freq file")
    }

    if (verbose) print(sprintf("Reading %s...", freqFileName))
    df = read.table(freqFileName, sep = ' ', stringsAsFactors = FALSE)

    markers = df[df[,1] == 'M',2]
    i = 1
    nMarker = 1
    for (markerName in markers) {
      i = i + 1
      als = c()
      freqs = c()
      while (i <= nrow(df) && df[i,1] == 'A') {
        als = c(als, df[i, 2])
        freqs = c(freqs, df[i, 3])
        i = i + 1
      }

      attr(ped$markerdata[[nMarker]], 'afreq') = freqs
      # TODO: assigning alleles causes genotypes to be changed kthxbye
      attr(ped$markerdata[[nMarker]], 'alleles') = als
      nMarker = nMarker + 1
    }
  }

  datFileName = paste(prefix, '.dat', sep = '')
  markernames = as.vector(lapply(ped$markerdata, function (m) { attr(m, 'name') }), mode = 'character')
  if ('dat' %in% what || (is.null(what) && length(markernames) > 0 && file.exists(datFileName) && file.size(datFileName) > 0)) {
    if (!file.exists(datFileName)) {
      stop2("Could not find requested .dat file")
    }

    if (verbose) print(sprintf("Reading %s...", datFileName))
    df = read.table(datFileName, sep = '', stringsAsFactors = FALSE)

    # TODO: use stopifnot2?
    stopifnot(df[,2] == markernames)
  }

  ped
}

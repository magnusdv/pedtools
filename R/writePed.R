#' Write a pedigree to file
#'
#' @param x A `ped` object
#' @param prefix A character string giving the prefix of the files. For
#'   instance, if `prefix = "myped"` and `what = c("ped", "map")`, the output
#'   files are "myped.ped" and "myped.map" in the current directory. Paths to
#'   other folder may be included, e.g. `prefix = "path-to-my-dir/myped"`.
#' @param what A subset of the character vector `c("ped", "map", "dat",
#'   "freq")`, indicating which files should be created. By default only the
#'   "ped" file is created. This option is ignored if `merlin = TRUE`.
#' @param famid A logical indicating if family ID should be included as the
#'   first column in the ped file. The family ID is taken from `famid(x)`. If
#'   `x` is a pedlist, the family IDs are taken from `names(x)`, or if this is
#'   NULL, the component-wise `famid()` values. Missing values are replaced by
#'   natural numbers. This option is ignored if `merlin = TRUE`.
#' @param header A logical indicating if column names should be included in the
#'   ped file. This option is ignored if `merlin = TRUE`.
#' @param merlin A logical. If TRUE, "ped", "map", "dat" and "freq" files are
#'   written in a format readable by the MERLIN software. In particular MERLIN
#'   requires non-numerical allele labels in the frequency file.
#' @param verbose A logical.
#'
#' @return A character vector with the file names.
#' @examples
#'
#' x = nuclearPed(1)
#' x = addMarker(x, "3" = "a/b", name = "m1")
#'
#' # Write to file
#' fn = writePed(x, prefix = tempfile("test"))
#'
#' # Read
#' y = readPed(fn)
#'
#' stopifnot(identical(x, y))
#'
#' @importFrom utils write.table
#' @export
writePed = function(x, prefix, what = "ped", famid = is.pedList(x),
                    header = TRUE, merlin = FALSE, verbose = TRUE) {

  if(merlin)
    return(writeMerlin(x, prefix = prefix, verbose = verbose))

  fnames = setNames(paste(prefix, what, sep = "."), what)

  if(is.pedList(x)) {
    pedmatr = do.call(rbind, lapply(x, as.data.frame.ped))

    if(famid) {
      famids = names(x) %||% unlist(lapply(x, famid))
      if(any(miss <- famids == "" | is.na(famids)))
        famids[miss] = seq_along(which(miss))
      famvec = rep(famids, pedsize(x))
      pedmatr = cbind(famid = famvec, pedmatr)
    }

    x = x[[1]] # for later use
  }
  else {
    pedmatr = as.data.frame(x)
    if(famid) {
      fam = if(x$FAMID == "") "1" else x$FAMID
      pedmatr = cbind(famid = fam, pedmatr)
    }
  }

  if ("ped" %in% what) {
    # TODO: This is slow; should use matrix objects.
    write.table(pedmatr, file = fnames[["ped"]], sep = "\t", col.names = header, row.names = FALSE, quote = FALSE)
    if(verbose) message("File written: ", fnames[["ped"]])

    # Faster, but excludes column names
    # write(t.default(pedmatr), fnames[["ped"]], ncolumns = ncol(pedmatr))
  }

  if ("freq" %in% what) {
    writeFreqDatabase(x, filename = fnames[["freq"]], format = "list")
    if(verbose) message("File written: ", fnames[["freq"]])
  }

  if (any(c("map", "dat") %in% what)) {
    mapmatr = getMap(x, na.action = 1, verbose = FALSE)
  }

  if ("map" %in% what) {
    write.table(mapmatr, fnames[["map"]], col.names = FALSE, row.names = FALSE, quote = FALSE)
    if(verbose) message("File written: ", fnames[["map"]])
  }

  if ("dat" %in% what) {
    datmatr = cbind("M", mapmatr$MARKER)
    write.table(datmatr, fnames[["dat"]], col.names = FALSE, row.names = FALSE, quote = FALSE)
    if(verbose) message("File written: ", fnames[["dat"]])
  }

  invisible(unname(fnames))
}



#' @importFrom utils write.table
writePed_merlin = function(x, prefix, verbose = TRUE) {

  what = c("ped", "map", "dat", "freq")
  fnames = setNames(paste(prefix, what, sep = "."), what)

  ### ped file
  if(is.pedList(x)) {
    pedmatr = do.call(rbind, lapply(x, as.matrix.ped, include.attrs = FALSE))
    pedmatr = cbind(rep.int(seq_along(x), pedsize(x)), pedmatr)
    x = x[[1]]
  } else {
    pedmatr = cbind(1, as.matrix(x, include.attrs = FALSE))
  }

  write(t.default(pedmatr), file = fnames[["ped"]], ncolumns = ncol(pedmatr))
  if(verbose) message("File written: ", fnames[["ped"]])

  ### map file
  mapmatr = getMap(x, na.action = 1, verbose = FALSE)
  write.table(mapmatr, file = fnames[["map"]], col.names = FALSE, row.names = FALSE, quote = FALSE)
  if(verbose) message("File written: ", fnames[["map"]])

  ### dat file
  datmatr = cbind("M", mapmatr$MARKER)
  write.table(datmatr, file = fnames[["dat"]], col.names = FALSE, row.names = FALSE, quote = FALSE)
  if(verbose) message("File written: ", fnames[["dat"]])

  ### freq file
  nalls = nAlleles(x)
  L = sum(nalls) + length(nalls)
  cum = cumsum(c(1, nalls + 1))
  length(cum) = length(nalls)  #remove last

  col1 = rep("A", L)
  col1[cum] = "M"

  col2 = character(L)
  col2[cum] = mapmatr$MARKER

  allalleles = unlist(lapply(nalls, seq_len)) # numerical allele names for merlin!
  col2[-cum] = allalleles

  col3 = character(L)
  allfreqs = unlist(lapply(x$MARKERS, afreq))
  col3[-cum] = format(allfreqs, scientifit = FALSE, digits = 6)

  freqmatr = cbind(col1, col2, col3)
  write.table(freqmatr, file = fnames[["freq"]], col.names = FALSE, row.names = FALSE, quote = FALSE)
  if(verbose) message("File written: ", fnames[["freq"]])

  invisible(unname(fnames))
}

#' @importFrom utils write.table
writeMerlin = function(x, prefix, verbose = TRUE) {

  what = c("ped", "map", "dat", "freq")
  fnames = setNames(paste(prefix, what, sep = "."), what)

  ### ped file
  if(is.pedList(x)) {
    pedmatr = do.call(rbind, lapply(x, as.matrix.ped, include.attrs = FALSE))
    pedmatr = cbind(rep.int(seq_along(x), pedsize(x)), pedmatr)
    x = x[[1]]
  } else {
    pedmatr = cbind(1L, as.matrix.ped(x, include.attrs = FALSE))
  }

  # writelines faster than write and write.table
  nr = nrow(pedmatr)
  lines = vapply(seq_len(nr), function(i) paste(pedmatr[i, ], collapse = " "), "")
  writeLines(lines, fnames[["ped"]])
  if(verbose) message("File written: ", fnames[["ped"]])

  ### map file
  mapmatr = getMap(x, merlin = TRUE, verbose = FALSE)
  write.table(mapmatr, file = fnames[["map"]], col.names = FALSE, row.names = FALSE, quote = FALSE)
  if(verbose) message("File written: ", fnames[["map"]])

  ### dat file
  markernames = mapmatr[,"MARKER"]
  datmatr = cbind("M", markernames)
  write.table(datmatr, file = fnames[["dat"]], col.names = FALSE, row.names = FALSE, quote = FALSE)
  if(verbose) message("File written: ", fnames[["dat"]])

  ### freq file
  nalls = nAlleles.ped(x)
  L = sum(nalls) + length(nalls)
  cum = cumsum(c(1, nalls + 1))
  length(cum) = length(nalls)  #remove last

  col1 = rep("A", L)
  col1[cum] = "M"

  col2 = character(L)
  col2[cum] = markernames

  allalleles = unlist(lapply(nalls, seq_len)) # numerical allele names for merlin!
  col2[-cum] = allalleles

  col3 = character(L)
  allfreqs = unlist(lapply(x$MARKERS, attr, "afreq"))
  col3[-cum] = format(allfreqs, scientifit = FALSE, digits = 6)

  freqmatr = cbind(col1, col2, col3)
  write.table(freqmatr, file = fnames[["freq"]], col.names = FALSE, row.names = FALSE, quote = FALSE)
  if(verbose) message("File written: ", fnames[["freq"]])

  invisible(unname(fnames))
}


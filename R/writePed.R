#' Write a pedigree to file
#'
#' @param x A `ped` object
#' @param prefix A character string giving the prefix of the files. For
#'   instance, if `prefix="myped"` and `what=c("ped", "map")`, the output files
#'   are "myped.ped" and "myped.map" in the current directory. Paths to other
#'   folder may be included, e.g. `prefix = "path-to-my-dir/myped"`. By default,
#'   the family identifier `famid(x)` is used. Note that if this is empty, the
#'   files become ".ped" and so on.
#' @param what A subset of the character vector c("ped", "map", "dat" and
#'   "freq"), indicating which files should be created. All files are written in
#'   MERLIN style (but see the `merlin` parameter below!) By default all files
#'   are created.
#' @param merlin A logical. If TRUE, the marker alleles are relabeled to
#'   1,2,..., making sure that the generated files are readable by MERLIN (which
#'   does not accept non-numerical allele labels in the frequency file.) If
#'   FALSE (the default) the allele labels are unchanged. In this case, `x`
#'   should be exactly reproducible from the files (see examples).
#' @param verbose A logical.
#'
#' @return A character vector with the file names.
#' @examples
#'
#' tmpdir = tempdir()
#' x = nuclearPed(1)
#' x = setMarkers(x, marker(x, '3'=1:2))
#' writePed(x, prefix = file.path(tmpdir, "myped"))
#'
#' @seealso [readPed()]
#'
#' @importFrom utils write.table
#' @export
writePed = function(x, prefix = famid(x), what = c("ped", "map", "dat", "freq"), merlin = FALSE, verbose=TRUE) {
  generated.files = character(0)

  if (merlin) {
    pedmatr = cbind(1, as.matrix(x, FALSE)) # 1 is FAMID
  }
  else { #TODO: avoid actual data.frame here. Slow! (Paramlink original ok...)
    famid = if(x$FAMID == "") 1 else x$FAMID
    pedmatr = cbind(famid=famid, as.data.frame(x))
  }

  if ("ped" %in% what) {
    pedname = paste(prefix, "ped", sep = ".")
    write(t(pedmatr), pedname, ncolumns = ncol(pedmatr))
    generated.files = c(generated.files, pedname)
    if(verbose) message("File written: ", pedname)
  }

  if (any(c("map", "dat", "freq") %in% what)) {
    mapmatr = getMap(x, na.action = 1, verbose = F)
    markerdata = x$markerdata
  }

  if ("map" %in% what) {
    mapname = paste(prefix, "map", sep = ".")
    write.table(mapmatr, mapname, col.names = F, row.names = F, quote = F)
    generated.files = c(generated.files, mapname)
    if(verbose) message("File written: ", mapname)
  }

  if ("dat" %in% what) {
    datname = paste(prefix, "dat", sep = ".")
    #datmatr = cbind(code = c("A", rep("M", nrow(mapmatr))), value = c("my_disease", mapmatr$MARKER))
    datmatr = cbind("M", mapmatr$MARKER)
    write.table(datmatr, datname, col.names = F, row.names = F, quote = F)
    generated.files = c(generated.files, datname)
    if(verbose) message("File written: ", datname)
  }

  if ("freq" %in% what) {
    freqname = paste(prefix, "freq", sep = ".")

    nalls = vapply(markerdata, nAlleles, 1)
    L = sum(nalls) + length(nalls)
    cum = cumsum(c(1, nalls + 1))
    length(cum) = length(nalls)  #remove last

    col1 = rep("A", L)
    col1[cum] = "M"

    col2 = character(L)
    col2[cum] = mapmatr$MARKER

    if (merlin)
      allalleles = unlist(lapply(nalls, seq_len))
    else
      allalleles = unlist(lapply(markerdata, alleles))

    col2[-cum] = allalleles

    col3 = character(L)
    allfreqs = unlist(lapply(markerdata, afreq))
    col3[-cum] = format(allfreqs, scientifit = F, digits = 6)

    freqmatr = cbind(col1, col2, col3)
    write.table(freqmatr, freqname, col.names = F, row.names = F, quote = F)
    generated.files = c(generated.files, freqname)
    if(verbose) message("File written: ", freqname)
  }

  invisible(generated.files)
}

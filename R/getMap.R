#' Tabulate marker positions
#'
#' Return a map of the markers attached to a pedigree.
#'
#' The `na.action` argument controls how missing values are dealt with:
#'
#' * `na.action` = 0: Return map unmodified
#'
#' * `na.action` = 1: Replace missing values with dummy values.
#'
#' * `na.action` = 2: Remove markers with missing data.
#'
#' In `setMap()`, the `map` argument should be a data frame (or file) with the
#' following columns in order: 1) chromosome, 2) marker name, 3) position in megabases.
#' Column names are ignored, as are any columns after the first three.
#'
#' @param x An object of class `ped` or a list of such.
#' @param markers A vector of names or indices referring to markers attached to
#'   `x`. By default, all markers are included.
#' @param na.action Either 0 (default), 1 or 2. (See Details.)
#' @param verbose A logical.
#' @param map Either a data frame or the path to a map file.
#' @param matchNames A logical; if TRUE, pre-existing marker names of `x` will
#'   be used to assign chromosome labels and positions from `map`.
#' @param ... Further arguments passed to `read.table()`.
#'
#' @return `getMap()` returns a data frame with columns `CHROM`, `MARKER` and
#'   `MB`.
#'
#'   `setMap()` returns `x` with modified marker attributes.
#'
#' @examples
#' x = singleton(1)
#' m1 = marker(x, chrom = 1, posMb = 10, name = "m1")
#' m2 = marker(x, chrom = 1, posMb = 11)
#' m3 = marker(x, chrom = 1)
#' x = setMarkers(x, list(m1, m2, m3))
#'
#' # Compare effect of `na.action`
#' getMap(x, na.action = 0)
#' getMap(x, na.action = 1)
#' getMap(x, na.action = 2)
#'
#' # Getting and setting map are inverses
#' y = setMap(x, getMap(x))
#' identical(x,y)
#'
#' @export
getMap = function(x, markers = NULL, na.action = 0, verbose = TRUE) {
  if(is.pedList(x)) {
    if(verbose) message("Input is a list of pedigrees; extracting map from first component")
    x = x[[1]]
  }
  if(!is.ped(x))
    stop2("Input is not a `ped` object or a list of such")

  if(!hasMarkers(x)) {
    if(verbose) message("No markers found")
    return(NULL)
  }

  if(is.null(markers))
    m = x$MARKERS
  else
    m = getMarkers(x, markers)

  chrom  = unlist(lapply(m, attr, "chrom"))
  marker = unlist(lapply(m, attr, "name"))
  mb     = unlist(lapply(m, attr, "posMb"))

  if (na.action == 1) {
    if (verbose)
      message("Warning: Missing map entries. Inserting dummy values.")

    naPos = is.na(mb) | is.na(chrom)

    # Autosomal unknowns: Put each on separate chromosome
    if(any(naPosAut <- naPos & !chrom %in% c("X", 23))) {
      chrom[naPosAut] = 100 + 1:sum(naPosAut)
      mb[naPosAut] = 0
    }

    # X unknowns: Separate by 400 cM
    if(any(naPosX <- naPos & chrom %in% c("X", 23))) {
      if(sum(naPosX) > 10) stop2("More than 10 markers on X with unknown position")
      mb[naPosX] = (1:sum(naPosX)) * 400 + if(all(naPos)) -400 else max(mb, na.rm = TRUE)
    }

    # NA names: make unique
    marker[is.na(marker)] = paste0("NA_", seq_len(sum(is.na(marker))))
  }
  else if(na.action == 2) {
    if(verbose)
      message("Warning: Missing map entries. Deleting markers with missing data.")
    miss = is.na(chrom) | is.na(marker) | is.na(mb)
    chrom = chrom[!miss]
    marker = marker[!miss]
    mb = mb[!miss]
  }

  data.frame(CHROM = chrom, MARKER = marker, MB = mb, stringsAsFactors = FALSE)
}



#' @rdname getMap
#' @importFrom utils read.table
#' @export
setMap = function(x, map, matchNames = NA, ...) {
  if(!is.ped(x) && !is.pedList(x))
    stop2("Input must be a `ped` object or a list of such")

  N = nMarkers(x)
  if(N == 0)
    stop2("The pedigree has no attached markers")

  if(is.character(map) && length(map) == 1)
    map = read.table(map, header = TRUE, as.is = TRUE, ...)

  if(!is.data.frame(map))
    stop2("`map` must be a data frame or file path")

  if(is.pedList(x)) {
    return(lapply(x, function(comp) setMap(comp, map = map, matchNames = matchNames)))
  }

  mapNames = map[[2]]
  xNames = name(x, 1:N)

  # Match names if either i) mismatch in number, or ii) names actually match in some order
  if(is.na(matchNames))
    matchNames = (nrow(map) != N) || (!all(is.na(mapNames)) && setequal(mapNames, xNames))

  if(matchNames) {
    mIdx = match(xNames, mapNames, nomatch = NA)
    mIdx = mIdx[!is.na(mIdx)]

    chrom(x, mIdx) = map[[1]][mIdx]
    posMb(x, mIdx) = map[[3]][mIdx]
  }
  else {
    if(nrow(map) != N)
      stop2("`map` incompatible with `x` (with `matchNames = F`)")
    chrom(x, 1:N) = map[[1]]
    posMb(x, 1:N) = map[[3]]
  }

  x
}

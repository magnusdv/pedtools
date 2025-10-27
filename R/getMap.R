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
#' following columns in order:
#'
#' 1) chromosome
#'
#' 2) marker name
#'
#' 3) position (Mb)
#'
#' Column names are ignored, as are any columns after the first three.
#'
#' @param x An object of class `ped` or a list of such.
#' @param markers A vector of names or indices referring to markers attached to
#'   `x`. By default, all markers are included.
#' @param na.action Either 0 (default), 1 or 2. (See Details.)
#' @param merlin A logical mostly for internal use: If TRUE the function returns
#'   a matrix instead of a data frame.
#' @param verbose A logical.
#' @param map Either a data frame, the path to a map file, or NULL (for removing
#'   map info). See Details regarding format.
#' @param matchNames A logical; if TRUE, pre-existing marker names of `x` will
#'   be used to assign chromosome labels and positions from `map`.
#' @param ... Further arguments passed to `read.table()`.
#'
#' @return `getMap()` returns a data frame with columns `CHROM`, `MARKER` and
#'   `MB`.
#'
#'   `setMap()` returns `x` with modified marker attributes.
#'
#'   `hasLinkedMarkers()` returns TRUE if two markers are located (with set
#'   position) on the same chromosome, and FALSE otherwise.
#'
#' @examples
#' x = singleton(1) |>
#'   addMarker(chrom = 1, posMb = 10, name = "m1") |>
#'   addMarker(chrom = 1, posMb = 11) |>
#'   addMarker(chrom = 1)
#'
#' # Compare effect of `na.action`
#' \donttest{
#' getMap(x, na.action = 0)
#' getMap(x, na.action = 1)
#' getMap(x, na.action = 2)
#' }
#'
#' # Getting and setting map are inverses
#' y = setMap(x, getMap(x))
#' stopifnot(identical(x,y))
#'
#' hasLinkedMarkers(x)
#'
#' @export
getMap = function(x, markers = NULL, na.action = 0, merlin = FALSE, verbose = TRUE) {
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

  if(merlin) {
    map = c(chrom, marker, mb)
    dim(map) = c(length(m), 3L)
    colnames(map) = c("CHROM", "MARKER", "MB")
    return(fixMerlinMap(map))
  }

  if(na.action == 1) {
    naPos = is.na(mb) | is.na(chrom)

    # Autosomal unknowns: Put each on separate chromosome
    if(any(naPosAut <- naPos & !chrom %in% c("X", 23))) {
      if(verbose)
        message("Warning: Missing map entries. Inserting dummy values.")

      chrom[naPosAut] = 100 + 1:sum(naPosAut)
      mb[naPosAut] = 0
    }

    # X unknowns: Separate by 400 cM
    if(any(naPosX <- naPos & chrom %in% c("X", 23))) {
      if(sum(naPosX) > 10) stop2("More than 10 markers on X with unknown position")
      if(verbose)
        message("Warning: Missing map entries. Inserting dummy values.")

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

  data.frame(CHROM = chrom, MARKER = marker, MB = mb)
}



# Internal merlin fix
fixMerlinMap = function(map) {
  if(!anyNA(map))
    return(map)

  chr = map[, "CHROM"]
  mar = map[, "MARKER"]
  mb = map[, "MB"]

  # NA names: make unique
  if(anyNA(mar))
    map[is.na(mar), "MARKER"] = paste0("_m", seq_len(sum(is.na(mar))))

  naPos = is.na(mb) | is.na(chr)
  if(!any(naPos))
    return(map)

  X = !is.na(chr) & (chr == "23" | chr == "X")
  if(!any(X)) Xchrom = FALSE
  else if(all(X)) Xchrom = TRUE
  else stop2("Mix of autosomal and X markers")

  # Autosomal: put on different chromosomes
  if(!Xchrom) {
    map[naPos, "CHROM"] = 100 + 1:sum(naPos)
    map[naPos, "MB"] = 0
  }
  else {  # X: Separate by 400 cM
    if(sum(naPos) > 10)
      stop2("More than 10 markers on X with unknown position")
    map[naPos, "MB"] = (1:sum(naPos)) * 400 + if(all(naPos)) -400 else max(mb, na.rm = TRUE)
  }

  map
}



#' @rdname getMap
#' @importFrom utils read.table
#' @export
setMap = function(x, map, matchNames = NA, ...) {

  N = nMarkers(x)
  if(N == 0)
    stop2("The pedigree has no attached markers")

  # Read in map if file path given
  if(is.character(map) && length(map) == 1)
    map = read.table(map, header = TRUE, as.is = TRUE, ...)

  if(is.pedList(x)) {
    return(lapply(x, function(comp) setMap(comp, map = map, matchNames = matchNames)))
  }

  ### Connected `ped` from here

  # If NULL map, set chrom and pos to NA and return
  if(is.null(map)) {
    return(x |> setChrom(chrom = NA) |> setPosition(posMb = NA))
  }

  if(!is.data.frame(map))
    stop2("`map` must be a data frame, a file path, or NULL")

  mapNames = map[[2]]
  xNames = name(x, 1:N)

  # Match names if either i) mismatch in number, or ii) names actually match in some order
  if(is.na(matchNames))
    matchNames = (nrow(map) != N) || (!any(is.na(mapNames)) && setequal(mapNames, xNames))

  if(matchNames) {
    mIdx = match(xNames, mapNames, nomatch = NA)
    mIdx = mIdx[!is.na(mIdx)]

    chr = map[[1]][mIdx]
    pos = map[[3]][mIdx]
    x = setChrom(x, 1:N, chrom = chr)
    x = setPosition(x, 1:N, posMb = pos)
  }
  else {
    if(nrow(map) != N)
      stop2("`map` incompatible with `x`. If the markers are named, set `matchNames = TRUE`")

    x = setChrom(x, 1:N, chrom = map[[1]])
    x = setMarkername(x, 1:N, name = map[[2]])
    x = setPosition(x, 1:N, posMb = map[[3]])
  }

  x
}



#' @export
#' @rdname getMap
hasLinkedMarkers = function(x) {
  map = getMap(x, na.action = 0, verbose = FALSE)
  if(is.null(map))
    return(FALSE)

  # With data on both chrom and position
  haspos = !is.na(map$CHROM) & !is.na(map$MB)

  # Return TRUE if two markers on the same chromosome
  anyDuplicated(map$CHROM[haspos]) > 0
}


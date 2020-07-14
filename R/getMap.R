
#' Tabulate marker positions
#'
#' @param x An object of class `ped`.
#' @param markers A numeric of indices.
#' @param pos Which unit should be used? Either "cm" (centiMorgan) or "mb" (megabytes).
#' @param na.action Either 0 (default), 1 or 2.
#' @param verbose A logical.
#'
#' @return A `data.frame`.
#' @export
#'
getMap = function(x, markers = seq_len(nMarkers(x)), pos = c("cm", "mb"), na.action = 0, verbose = TRUE) {
  # TODO review this function
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(!hasMarkers(x)) return(NULL)
  m = getMarkers(x, markers)
  chrom = unlist(lapply(m, attr, "chrom"))
  marker = unlist(lapply(m, attr, "name"))
  pos = switch(match.arg(pos),
               cm = unlist(lapply(m, attr, "posCm")),
               mb = unlist(lapply(m, attr, "posMb")))
  map = data.frame(CHR = chrom, MARKER = marker, POS = pos, stringsAsFactors = FALSE)
  if (na.action > 0) {
    na_pos = (is.na(chrom) | is.na(pos))
    na_name = is.na(marker)
    map$MARKER[na_name] = paste0("M", markers[na_name])
  }
  if (na.action == 1 && all(na_pos)) {
    if (verbose) message("Warning: No map info given. Creating dummy map.")
    map$CHR = rep_len(1, nrow(map))
    map$POS = seq_len(nrow(map))
  }
  if(na.action == 2 && any(na_pos)) {
    if(verbose)
      message('Warning: Deleting ', sum(na_pos),
              ' markers with missing map coordinates.')
    map = map[!na_pos, , drop = FALSE]
  }
  map
}



#' @importFrom utils read.table
#' @export
setMap = function(x, map, matchNames = NA, ...) {
  if(!is.ped(x) || is.pedList(x))
    stop2("Input must be a `ped` object or a list of such")

  N = nMarkers(x)
  if(N == 0)
    stop2("The pedigree has no attached markers")

  if(is.character(map) && length(map) == 1)
    map = read.table(map, header = TRUE, as.is = TRUE, ...)

  if(!is.data.frame(map))
    stop2("`map` must be a data frame or file path")

  mapNames = map[[2]]
  xNames = name(x, 1:N)

  # Match names if either i) mismatch in number, or ii) names actually match in some order
  if(is.na(matchNames))
    matchNames = (nrow(map) != N) || setequal(mapNames, xNames)

  if(matchNames) {
    mIdx = match(xNames, mapNames, nomatch = NA)
    mIdx = mIdx[!is.na(mIdx)]

    chrom(x, mIdx) = map[[1]][mIdx]
    posMb(x, mIdx) = map[[3]][mIdx]
  }
  else {
    if(nrow(map) != N)
      stop2("`map` incompatible with `x` (with `matchNames = F`)")
    chrom(x, N) = map[[1]]
    name(x, N) = map[[3]]
    posMb(x, mIdx) = map[[3]]
  }

  x
}


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
setMap = function(x, map) {
  if(is.character(map))
    map = read.table(map, header = TRUE)

  chroms = map[[1]]
  mnames = map[[2]]
  chrom(x, mnames) = chroms

  MBcol = grep("mb", names(map), value = TRUE, ignore.case = TRUE)
  CMcol = grep("cm", names(map), value = TRUE, ignore.case = TRUE)

  if(length(MBcol))
    posMb(x, mnames) = map[[MBcol[1]]]
  if(length(CMcol))
    posCm(x, mnames) = map[[CMcol[1]]]

  x
}

#' Attach identical markers evenly distributed across the genome
#'
#' @param x A ped object.
#' @param n The total number of markers. Either this or `dist` must be NULL.
#' @param dist A positive number; the distance (in megabases) between markers.
#' @param chromLen A numeric vector indicating chromosome lengths (in Mb). By
#'   default, the lengths of the human chromosomes 1-22 are used, as returned by
#'   `sapply(ibdsim2::loadMap("decode"), ibdsim2::physRange)`.
#' @param alleles,afreq Passed onto [marker()].
#' @param prefix A string used as prefix for marker names. The default (NULL)
#'   gives no names.
#'
#' @return A ped object similar to `x`.
#'
#' @examples
#' x = distributeMarkers(nuclearPed(), n = 10, prefix = "M")
#' getMap(x)
#'
#' @export
distributeMarkers = function(x, n = NULL, dist = NULL, chromLen = NULL,
                             alleles = 1:2, afreq = NULL, prefix = NULL) {
  if(is.null(chromLen))
    chromLen = c(246.98258, 241.01465, 197.08100, 188.84446, 180.22043, 169.42632, 158.18255,
                            143.83232, 136.96897, 132.52836, 133.86230, 132.11300, 94.68256, 85.03331, 77.41597,
                            89.03333, 81.90883, 79.11141, 57.31905, 63.19036, 32.01867, 33.16128)
  # Total genome length
  L = sum(chromLen)

  if(!is.null(n) && !is.null(dist))
    stop2("Exactly one of `n` and `dist` must be given")

  # Compute positions if only the number of markers is given
  if(!is.null(n)) {
    pos0 = seq(0, L, length.out = n)
    starts0 = cumsum(c(0, chromLen))
    posList = split(pos0, cut(pos0, starts0, labels = FALSE, include.lowest = TRUE))
    chr = as.integer(names(posList))
    starts = starts0[chr]
    pos = unlist(unname(posList)) - rep(starts, lengths(posList))
  }
  else {
    posList = lapply(chromLen, function(len) seq(from = 0, to = len, by = dist))
    pos = unlist(posList, use.names = FALSE)
    chr = rep(seq_along(posList), lengths(posList))
  }


  nms = if(!is.null(prefix)) paste0(prefix, seq_along(chr)) else NA_character_  # note: NA accepts subsetting

  m = marker(x, alleles = afreq, afreq = afreq)
  mlist = lapply(seq_along(pos), function(i) {
    mi = m
    attr(mi, "chrom") = chr[i]
    attr(mi, "posMb") = pos[i]
    attr(mi, "name") = nms[i]
    mi
  })

  # Speeds up setMarkers()
  class(mlist) = "markerList"

  # Attach and return
  setMarkers(x, mlist, checkCons = FALSE)
}



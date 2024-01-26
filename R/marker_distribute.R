#' Distribute markers evenly along a set of chromosomes
#'
#' Create and attach identical (empty) marker objects, distributed along a set
#' of chromosomes.
#'
#' Note: When using the `dist` parameter, the function treats each chromosome
#' separately, places one marker at the start and then every `dist` megabases.
#' (See Examples.)
#'
#' @param x A ped object.
#' @param n The total number of markers. Either this or `dist` must be NULL.
#' @param dist A positive number; the distance (in megabases) between markers.
#' @param chromLen A numeric vector indicating chromosome lengths (in Mb). By
#'   default, the lengths of the human chromosomes 1-22 are used, as returned by
#'   `sapply(ibdsim2::loadMap("decode"), ibdsim2::physRange)`.
#' @param alleles,afreq Passed onto [marker()].
#' @param prefix A string used as prefix for marker names. Default: "M".
#'
#' @return A copy of `x` with the indicated markers attached.
#'
#' @examples
#' x = distributeMarkers(nuclearPed(), n = 10)
#' getMap(x)
#'
#' y = distributeMarkers(nuclearPed(), dist = 100)
#' getMap(y)
#' @export
distributeMarkers = function(x, n = NULL, dist = NULL, chromLen = NULL,
                             alleles = 1:2, afreq = NULL, prefix = "M") {

  if(!length(chromLen)) {
    chromLen = c(246.98258, 241.01465, 197.08100, 188.84446, 180.22043, 169.42632, 158.18255,
                 143.83232, 136.96897, 132.52836, 133.86230, 132.11300, 94.68256, 85.03331, 77.41597,
                 89.03333, 81.90883, 79.11141, 57.31905, 63.19036, 32.01867, 33.16128)
  }
  chromNames = names(chromLen) %||% as.character(seq_along(chromLen))

  if(any(chromNames == ""))
    stop2("Irregular chromosome names")
  if(dup <- anyDuplicated(chromNames))
    stop2("Duplicated chromosome name: ", chromNames[dup])

  # Total genome length
  L = sum(chromLen)

  if(!is.null(n) && !is.null(dist))
    stop2("Exactly one of `n` and `dist` must be given")

  # Compute positions if only the number of markers is given
  if(!is.null(n)) {

    if(length(n) != 1 || !is.numeric(n) || n <= 0 || round(n) != n)
      stop2("`n` must be a positive integer (or NULL): ", n)

    pos0 = seq(0, L, length.out = n)
    starts0 = cumsum(c(0, chromLen))
    posList = split(pos0, cut(pos0, starts0, labels = FALSE, include.lowest = TRUE))
    chrnum = as.integer(rep(names(posList), lengths(posList)))
    starts = starts0[chrnum]
    pos = unlist(posList, use.names = FALSE) - starts
  }
  else {

    if(length(dist) != 1 || !is.numeric(dist) || dist <= 0)
      stop2("`dist` must be a positive number (or NULL): ", dist)

    posList = lapply(chromLen, function(len) seq(from = 0, to = len, by = dist))
    pos = unlist(posList, use.names = FALSE)
    chrnum = rep(seq_along(posList), lengths(posList))
  }

  # Chromosome names
  chr = chromNames[chrnum]

  # Marker names
  nms = if(!is.null(prefix)) paste0(prefix, seq_along(chr)) else NA_character_  # note: NA accepts subsetting

  # Alleles/freqs
  if(!is.null(names(afreq)))
    alleles = NULL

  m = marker(x, alleles = alleles, afreq = afreq)
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



#' Attach SNP loci to a pedigree
#'
#' Create and attach a list of empty SNP markers with specified position and
#' allele frequencies.
#'
#' The data frame `snpData` should contain the following columns, in order:
#'
#' * `CHROM`: Chromosome (character)
#'
#' * `MARKER`: Marker name (character)
#'
#' * `MB`: Physical position in megabases (numeric)
#'
#' * `A1`: First allele (single-letter character)
#'
#' * `A2`: Second allele (single-letter character)
#'
#' * `FREQ1`: Allele frequency of `A1` (number in `[0,1]`)
#'
#' The actual column names do not matter.
#'
#' Each column must be of the stated type, or coercible to it. (For example,
#' `CHROM`, `A1` and `A2` may be given as numbers, but will be internally
#' converted to characters.)
#'
#' @param x A `ped` object.
#' @param snpData A data frame with 6 columns. See Details.
#'
#' @return A copy of `x` with the indicated SNP markers attached.
#'
#' @examples
#'
#' snps = data.frame(
#'   CHROM  = 1:2,
#'   MARKER = c("M1", "M2"),
#'   MB     = c(1.23, 2.34),
#'   A1     = c("A", "G"),
#'   A2     = c("C", "C"),
#'   FREQ1  = c(0.7, 0.12))
#'
#' x = setSNPs(nuclearPed(), snpData = snps)
#'
#' # Inspect the results:
#' getMap(x)
#' getFreqDatabase(x)
#'
#' @importFrom utils head
#' @export
setSNPs = function(x, snpData) {

  if(is.pedList(x))
    return(lapply(x, function(y) setSNPs(y, snpData)))
  if(!is.ped(x))
    stop2("Argument `x` is not a ped object: ", class(x))

  if(!is.data.frame(snpData) || ncol(snpData) != 6)
    stop2("`snpData` must be a data frame with 6 columns. See ?setSNPs")

  names(snpData) = toupper(names(snpData))

  # If all names present, but possibly wrong order: sort
  nms = c("CHROM", "MARKER", "MB", "A1", "A2", "FREQ1")
  if(setequal(nms, names(snpData)))
    snpData = snpData[nms]

  # Columns
  CHROM = as.character(snpData[[1]])
  MARKER = as.character(snpData[[2]])
  MB = as.numeric(snpData[[3]])
  A1 = as.character(snpData[[4]])
  A2 = as.character(snpData[[5]])
  FREQ1 = as.numeric(snpData[[6]])

  if(any(bad <- (FREQ1 > 1 | FREQ1 < 0)))
    stop2("Illegal allele frequency: ", head(FREQ1[bad]))
  als = c(A1, A2)
  legal = c(LETTERS, letters, 1:9)
  if(!all(als %in% legal))
    stop2("Illegal allele label (should be single letters/digits): ", head(setdiff(als, legal)))

  # Sort allele pairs
  swap = A1 > A2
  if(any(swap)) {
    A1tmp = A1; A1[swap] = A2[swap]; A2[swap] = A1tmp[swap]
    FREQ1[swap] = 1 - FREQ1[swap]
  }

  # Same for all markers
  labs = labels(x)
  sex = getSex(x)
  amat = matrix(0L, ncol = 2, nrow = pedsize(x))

  # Construct markers
  mlist = lapply(seq_along(MARKER), function(i)
    newMarker(amat, alleles = c(A1[i], A2[i]), afreq = c(FREQ1[i], 1 - FREQ1[i]),
              name = MARKER[i], chrom = CHROM[i], posMb = MB[i], pedmembers = labs, sex = sex)
  )

  class(mlist) = "markerList"
  setMarkers(x, mlist, checkCons = FALSE)
}



.setSNPfreqs = function(x, freq1) {
  if(!is.numeric(freq1))
    stop2("Argument `freq1` must be numeric, not ", class(freq1)[1])

  if(any(bad <- (freq1 < 0 | freq1 > 1)))
    stop2("Element of `freq1` outside interval [0,1]: ", freq1[bad])

  nm = nMarkers(x)
  if(length(freq1) == 1)
    freq1 = rep(freq1, nm)
  else if(length(freq1) != nm)
    stop2("Length of `freq1` must equal the number of markers (or 1)")

  for(idx in seq_len(nMarkers(x))) {
    m = x$MARKERS[[idx]]

    # Check diallelic
    if(length(attr(m, "alleles")) != 2) {
      mname = attr(m, 'name')
      if(is.na(mname))
        mname = sprintf("<%d>", idx)
      stop2("Marker does not have exactly 2 alleles: ", mname)
    }

    fr = as.numeric(freq1[idx])
    attr(m, "afreq") = c(fr, 1 - fr)
    x$MARKERS[[idx]] = m

    if(allowsMutations(m))
      x = setMutmod(x, markers = idx, update = TRUE)
  }

  x
}

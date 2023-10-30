#' Mask and unmask pedigree datasets
#'
#' The `maskPed()` function replaces the individual IDs, marker names and allele
#' labels with generic sequences like 1, 2, ... For markers with stepwise
#' mutation models, the allelic ladder is simply translated to start at 1,
#' thereby preserving the intra-allelic differences.
#'
#' It should be noted that when the masking procedure is applied to a dataset
#' using publicly available frequency databases, the result cannot be considered
#' to be fully anonymised. (In theory, one could deduce the original marker
#' names and alleles from the frequencies.)
#'
#' @param x A `ped` object or a list of such.
#' @param seed An optional seed for the random number generator.
#' @param keys A list with entries `ids`, `markers`, `alleles`.
#'
#' @return An object similar to `x` but with replaced ID labels, marker names
#'   and allele labels.
#'
#' @examples
#' x = nuclearPed(father = "fa", mother = "mo", children = "ch") |>
#'   addMarker(name = "myMarker", ch = "b/c", afreq = c(a=0.2, b=0.3, c=0.5)) |>
#'   setMutmod(model = "proportional", rate = 0.01)
#'
#' # Mask
#' y = maskPed(x, seed = 1729)
#'
#' # Unmask
#' z = unmaskPed(y$maskedPed, y$keys)
#' stopifnot(identical(x, z))
#'
#' # With stepwise model
#' x2 = x |>
#'   addMarker(name = "mySTR", ch = "7.2/8.2",
#'             alleles = c("7", "7.2", "8", "8.2")) |>
#'   setMutmod(marker = 2, model = "stepwise", rate = 0.1, rate2 = 1e-6,
#'             range = 0.1)
#'
#' y2 = maskPed(x2, seed = 1729)
#'
#' z2 = unmaskPed(y2$maskedPed, y2$keys)
#'
#' stopifnot(identical(x2, z2))
#'
#' # Check likelihoods with pedprobr:
#' # stopifnot(setequal(likelihood(x2), likelihood(y2$maskedPed)))
#'
#' @importFrom stats setNames
#' @export
maskPed = function(x, seed = NULL) {

  if(!is.null(seed))
    set.seed(seed)

  # Individuals 1,2,...
  ids = unlist(labels(x), use.names = FALSE)
  newids = as.character(seq_along(ids))
  names(newids) = ids
  y = relabel(x, newids)

  nm = nMarkers(y)
  if(nm == 0)
    return(list(maskedPed = y, keys = list(ids = newids)))

  mnames = name(y)
  if(anyNA(mnames))
    stop2("The masking procedure requires all markers to be named")

  # Create new allele labels
  newAls = lapply(1:nm, function(i) {
    als = alleles(y, marker = i)
    alsNum = suppressWarnings(as.numeric(als))
    if(!any(is.na(alsNum)))
      nw = round(alsNum - min(alsNum) + 1, 1)
    else
      nw = sample.int(length(als))
    setNames(as.character(nw), als)
  })
  names(newAls) = name(y)

  # Apply new allele labels
  for(i in 1:nm)
    y = setAlleleLabels(y, marker = i, alleles = newAls[[i]])

  # Shuffle markers
  shuffledNames = sample(mnames)
  y = selectMarkers(y, markers = shuffledNames)

  # Marker names M1, M2, ...
  newnames = setNames(paste0("M", 1:nm), shuffledNames)
  y = setMarkername(y, name = newnames)

  # Sort `newnames` in original order
  newnames = newnames[mnames]

  # Return pedigree and keys for unmasking
  list(maskedPed = y, keys = list(ids = newids, alleles = newAls, markernames = newnames))
}

#' @rdname maskPed
#' @export
unmaskPed = function(x, keys) {

  # Restore individuals IDs
  y = relabel(x, .flipNames(keys$ids))

  # Restore marker names
  y = setMarkername(y, name = .flipNames(keys$markernames))

  # Restore marker order
  y = selectMarkers(y, names(keys$markernames))

  # Restore allele labels
  for(m in names(keys$alleles))
    y = setAlleleLabels(y, marker = m, alleles = .flipNames(keys$alleles[[m]]))

  y
}

# Swap elements and names of a vector
.flipNames = function(v) {
  setNames(names(v), v)
}

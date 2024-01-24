#' Mask and unmask pedigree datasets
#'
#' The `maskPed()` function replaces the individual IDs, marker names and allele
#' names with generic labels, and randomly changes their internal order. For
#' markers with stepwise mutation models, the allelic ladder is simply
#' translated to start at 1, thereby preserving the intra-allelic differences.
#'
#' Note that in order to preserve likelihoods, the allele frequencies are not
#' modified. Thus, if the data uses a publicly available frequency databases,
#' the result cannot be considered to be fully anonymised, since one could (at
#' least in theory) deduce the original marker names and alleles from the
#' frequencies.)
#'
#' @param x A `ped` object or a list of such.
#' @param ids (Optional) A named character with the new IDs, written as `c(old =
#'   new, ...)`. By default: 1, 2, ... .
#' @param markerNames (Optional) A named character with the new marker names
#'   (and order), written as `c(old = new, ...)`. By default: M1, M2, ... .
#' @param markerShuffle A logical: Randomly reorder the markers? (Default: TRUE)
#' @param alleleLabels (Optional) A list of character vectors. The list names
#'   should be the original marker names. Each vector gives the new allele
#'   labels, as `c(old = new, ...)`. By default, each marker gets alleles 1, 2,
#'   ... .
#' @param alleleShuffle A logical: Randomly reorder the alleles? (Default: TRUE)
#' @param seed An optional seed for the random number generator.
#' @param keys A list with entries `ids`, `markerNames`, `alleleLabels`.
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
#' z = unmaskPed(y$maskedPed, keys = y$keys)
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
#' z2 = unmaskPed(y2$maskedPed, keys = y2$keys)
#'
#' stopifnot(identical(x2, z2))
#'
#' # Check likelihoods with pedprobr:
#' # stopifnot(setequal(likelihood(x2), likelihood(y2$maskedPed)))
#'
#' @importFrom stats setNames
#' @export
maskPed = function(x, ids = NULL, markerNames = NULL, markerShuffle = TRUE,
                   alleleLabels = NULL, alleleShuffle = TRUE, seed = NULL) {

  if(!is.null(seed))
    set.seed(seed)

  # Individuals (default: 1,2,...)
  if(is.null(ids)) {
    oldids = unlist(labels(x), use.names = FALSE)
    ids = as.character(seq_along(oldids))
    names(ids) = oldids
  }
  y = relabel(x, ids)

  nm = nMarkers(y)
  if(nm == 0)
    return(list(maskedPed = y, keys = list(ids = ids)))

  oldnames = name(y)
  if(anyNA(oldnames))
    stop2("The masking procedure requires all markers to be named")

  # Allele labels
  if(is.null(alleleLabels)) {
    alleleLabels = lapply(1:nm, function(i) {
      alsOld = alleles(y, marker = i)
      alsNum = suppressWarnings(as.numeric(alsOld))
      if(!any(is.na(alsNum)) && .hasStepwiseModel(y, i))
        als = round(alsNum - round(min(alsNum)) + 1, 1)
      else if(alleleShuffle)
        als = sample.int(length(alsOld))
      else
        als = seq_along(alsOld)
      setNames(as.character(als), alsOld)
    })
    names(alleleLabels) = oldnames
  }

  # Apply new allele labels
  # NB: Do this before shuffling (much faster to supply idx than name)
  for(i in seq_along(oldnames))
    y = setAlleleLabels(y, marker = i, alleles = alleleLabels[[i]])

  # Marker order
  if(!is.null(markerNames))
    mOrder = names(markerNames)
  else if(markerShuffle)
    mOrder =  sample(oldnames)
  else
    mOrder = oldnames

  # Rename markers (default: M1, M2, ...)
  if(is.null(markerNames))
    markerNames = setNames(paste0("M", 1:nm), mOrder)

  y = selectMarkers(y, markers = mOrder)
  y = setMarkername(y, name = markerNames)

  # For key: sort markerNames in original order
  markerNames = markerNames[oldnames]

  # Return pedigree and keys for unmasking
  list(maskedPed = y, keys = list(ids = ids, alleleLabels = alleleLabels,
                                  markerNames = markerNames))
}

#' @rdname maskPed
#' @export
unmaskPed = function(x, keys) {

  # Restore individuals IDs
  y = relabel(x, .flipNames(keys$ids))

  # Restore marker names
  y = setMarkername(y, name = .flipNames(keys$markerNames))

  # Restore marker order
  y = selectMarkers(y, names(keys$markerNames))

  # Restore allele labels
  for(m in names(keys$alleleLabels))
    y = setAlleleLabels(y, marker = m, alleles = .flipNames(keys$alleleLabels[[m]]))

  y
}

# Swap elements and names of a vector
.flipNames = function(v) {
  setNames(names(v), v)
}

.hasStepwiseModel = function(x, marker) {
  mut = mutmod(x, marker = marker)
  if(is.null(mut))
    return(FALSE)
  params = pedmut::getParams(mut, format = 1)
  any(params$model %in% c("stepwise", "onestep"))
}

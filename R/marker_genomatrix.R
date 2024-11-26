#' Genotype matrix
#'
#' Extract the genotypes of specified individuals and markers from a pedigree
#' object, and return them as a character matrix.
#'
#' @param x A `ped` object or a list of such.
#' @param ids A vector of ID labels, or a function operating on `x`, e.g.,
#'   [typedMembers()]. By default (`ids = NULL`) all individuals are included,
#'   also non-genotyped ones.
#' @param markers A vector of indices or names of markers attaches to `x`. If
#'   NULL (default) all markers are included.
#' @param sep A single string to be used as allele separator in marker
#'   genotypes.
#' @param missing A single string to be used for missing alleles.
#'
#' @return
#'
#' `getGenotypes()` returns a character matrix with `length(ids)` rows and
#' `length(markers)` columns.
#'
#' @seealso [getAlleles()]
#'
#' @examples
#' x = nuclearPed() |>
#'   addMarker(`2` = "1/2", name = "m1") |>
#'   addMarker(`3` = "a/a", name = "m2")
#'
#' getGenotypes(x)
#'
#' ### A list of pedigrees
#'
#' s = transferMarkers(x, singleton("s"))
#' peds = list(x, s)
#'
#' getGenotypes(peds)
#'
#' # Using a function to select individuals
#' getGenotypes(x, ids = typedMembers)
#'
#' @export
getGenotypes = function(x, ids = NULL, markers = NULL, sep = "/", missing = "-") {
  discon = !is.ped(x)
  if(discon && !is.pedList(x)) {
    stop2("The first argument must be a `ped` object or a list of such")
  }

  # Check `ids` argument
  if(!is.null(ids)) {

    if(is.function(ids))
      ids = ids(x)
    else {
      ids = as.character(ids)
      if(!all(ids %in% labels(x)))
        stop2("Unknown ID label: ", setdiff(ids, labels(x)))
      if(dup <- anyDuplicated.default(ids))
        stop2("Duplicated element of argument `ids`: ", dup)
    }
  }

  if(discon) {

    # Check equality of marker counts and names
    name(x, seq_markers(x))

    # Extract genotypes from each component
    compList = lapply(x, function(comp) {
      ids_comp = if(!is.null(ids)) .myintersect(ids, labels(comp))
      getGenotypes(comp, ids = ids_comp, markers = markers)
    })

    # Bind into single matrix
    res = do.call(rbind, compList)

    # Sort rows according to input `ids`
    if(!is.null(ids))
      res = res[ids, , drop = FALSE]

    return(res)
  }

  # Main body: x is now is single `ped` object
  ids = ids %||% x$ID

  # Ensure markers is integer
  markers = if(!is.null(markers)) whichMarkers(x, markers) else seq_markers(x)

  # Headers of genotype columns: name if present, otherwise <idx>
  mNames = name(x, markers)
  if(any(miss <- is.na(mNames)))
    mNames[miss] = sprintf("<%d>", markers[miss])

  # If no `ids` or no `markers`, return empty
  if(length(ids) == 0 || length(markers) == 0) {
    empt = matrix(character(0), nrow = length(ids), ncol = length(markers),
                  dimnames = list(ids, mNames))
    return(empt)
  }

  # Genotype matrix with all individuals
  mlist = getMarkers(x, markers)
  g = do.call(cbind, lapply(mlist, function(m) format(m, sep = sep, missing = missing)))

  # Set dimnames
  rownames(g) = labels(x)
  colnames(g) = mNames

  # Return subset
  g[ids, , drop = FALSE]  # NB: ids must be character here
}

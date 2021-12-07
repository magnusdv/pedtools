#' Genotype matrix
#'
#' Extract the genotypes of multiple individuals/markers in form of a matrix.
#'
#' @param x A `ped` object or a list of such
#' @param ids A vector of ID labels. If NULL (default) all individuals are
#'   included.
#' @param markers A vector of indices or names of markers attaches to `x`. If
#'   NULL (default) all markers are included.
#' @param sep A single string to be used as allele separator in marker
#'   genotypes.
#' @param missing A single string to be used for missing alleles.
#'
#' @return
#'
#' `getGenotypes()` returns a character matrix with `length(ids)` rows
#'   and `length(markers)` columns.
#'
#' @seealso [getAlleles()]
#'
#' @examples
#' x = nuclearPed(1)
#' m1 = marker(x, `2` = "1/2", alleles = 1:2, name = "m1")
#' m2 = marker(x, `3` = "2/2", alleles = 1:2, name = "m2")
#' x = setMarkers(x, list(m1, m2))
#'
#' getGenotypes(x)
#'
#'
#' ### A list of pedigrees
#'
#' s = transferMarkers(x, singleton("s"))
#' peds = list(x, s)
#'
#' getGenotypes(peds)
#'
#' @export
getGenotypes = function(x, ids = NULL, markers = NULL, sep = "/", missing = "-") {
  if(!is.ped(x) && !is.pedList(x))
    stop2("The first argument must be a `ped` object or a list of such")

  if(dup <- anyDuplicated(ids)) {
    stop2("Duplicated element of argument `ids`: ", dup)
  }

  # Check `ids` argument
  if(!is.null(ids)) {
    ids = as.character(ids)
    if(!all(ids %in% unlist(labels(x))))
      stop2("Unknown ID label: ", setdiff(ids, unlist(labels(x))))
  }

  if(is.pedList(x)) {

    # Check equality of marker counts and names
    name(x, seq_markers(x))

    # Extract genotypes from each component
    compList = lapply(x, function(comp) {
      ids_comp = if(!is.null(ids)) intersect(ids, labels(comp))
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
  ids = ids %||% labels(x)

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

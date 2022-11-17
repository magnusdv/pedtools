#' Allele matrix manipulation
#'
#' Functions for getting and setting the genotypes of multiple
#' individuals/markers simultaneously
#'
#' If the `alleles` argument of `setAlleles()` is not a matrix, it is recycled
#' (if necessary), and converted into a matrix of the correct dimensions. For
#' example, setting `alleles = 0` gives a simple way of removing the genotypes
#' of some or all individuals (while keeping the markers attached).
#'
#' @param x A `ped` object or a list of such
#' @param ids A vector of ID labels. If NULL (default) all individuals are
#'   included.
#' @param markers A vector of indices or names of markers attaches to `x`. If
#'   NULL (default) all markers are included.
#' @param alleles A character of the same format and dimensions as the output of
#'   `getAlleles(x, ids, markers)`, or an object which can be converted by
#'   `as.matrix()` into such a matrix. See Details.
#'
#' @return `getAlleles()` returns a character matrix with `length(ids)` rows and
#'   `2 * length(markers)` columns. The ID labels of `x` are used as rownames,
#'   while the columns are named `<m1>.1`, `<m1>.2`, ... where `<m1>` is the
#'   name of the first marker, a.s.o.
#'
#'   `setAlleles()` returns a `ped` object identical to `x`, except for the
#'   modified alleles. In particular, all locus attributes are unchanged.
#'
#' @seealso [transferMarkers()]
#'
#' @examples
#' # Setup: Pedigree with two markers
#' x = nuclearPed(1)
#' x = addMarker(x, `2` = "1/2", alleles = 1:2, name = "m1")
#' x = addMarker(x, `3` = "2/2", alleles = 1:2, name = "m2")
#'
#' # Extract allele matrix
#' mat1 = getAlleles(x)
#' mat2 = getAlleles(x, ids = 2:3, markers = "m2")
#' stopifnot(identical(mat1[2:3, 3:4], mat2))
#'
#' # Remove all genotypes
#' y = setAlleles(x, alleles = 0)
#' y
#'
#' # Setting a single genotype
#' z = setAlleles(y, ids = "1", marker = "m2", alleles = 1:2)
#'
#' # Alternative: In-place modification with `genotype()`
#' genotype(y, id = "1", marker = "m2") = "1/2"
#' stopifnot(identical(y,z))
#'
#'
#' ### Manipulation of pedlist objects
#' s = transferMarkers(x, singleton("s"))
#' peds = list(x, s)
#'
#' getAlleles(peds)
#'
#' setAlleles(peds, ids = "s", marker = "m1", alleles = 1:2)
#'
#' @export
getAlleles = function(x, ids = NULL, markers = NULL) {
  if(!is.ped(x) && !is.pedList(x))
    stop2("The first argument must be a `ped` object or a list of such")

  if(dup <- anyDuplicated(ids)) {
    stop2("Duplicated element of argument `ids`: ", dup)
  }

  if(is.pedList(x)) {

    # Check that all `ids` are known
    if(!is.null(ids) && !all(ids %in% unlist(labels(x))))
      stop2("Unknown ID label: ", setdiff(ids, unlist(labels(x))))

    # Check equality of marker counts and names
    name(x)

    # Extract alleles from each component
    amList = lapply(x, function(comp) {
      ids_comp = if(!is.null(ids)) intersect(ids, labels(comp))
      getAlleles(comp, ids = ids_comp, markers = markers)
    })

    # Bind into single matrix
    res = do.call(rbind, amList)

    # Sort rows according to input `ids`
    if(!is.null(ids))
      res = res[as.character(ids), , drop = FALSE]

    return(res)
  }

  # Main body: x is now is single `ped` object
  ids = ids %||% labels(x)
  markers = markers %||% seq_markers(x)

  # If no `ids` or no `markers`, return empty (but properly formed) matrix
  if(length(ids) == 0 || length(markers) == 0) {
    am = matrix(character(0), nrow = length(ids), ncol = 2*length(markers),
                dimnames = list(ids, sprintf("%s.%d", rep(markers, each = 2), 1:2)))
    return(am)
  }

  mlist = getMarkers(x, markers = markers)

  # Convert to character matrix
  am = markerlist2allelematrix(mlist)

  # Subset and set rownames
  amSubset = am[internalID(x, ids), , drop = FALSE]
  rownames(amSubset) = ids

  # Return
  amSubset
}

#' @rdname getAlleles
#' @export
setAlleles = function(x, ids = NULL, markers = NULL, alleles) {
  if(!is.ped(x) && !is.pedList(x))
    stop2("The first argument must be a `ped` object or a list of such")

  completeAlleleMatrix = getAlleles(x, markers = markers)

  if(is.null(ids))
    ids = rownames(completeAlleleMatrix)
  else {
    ids = as.character(ids)
    if(!all(ids %in% rownames(completeAlleleMatrix)))
      stop2("Unknown ID label: ", setdiff(ids, rownames(completeAlleleMatrix)))
  }

  oldAlleles = completeAlleleMatrix[ids, , drop = FALSE]
  if(is.null(oldAlleles))
    return(x)

  # Fix `alleles` if not matrix
  if(is.data.frame(alleles))
    alleles = as.matrix(alleles)

  if(length(alleles) == 1)
    alleles = rep(alleles, length(oldAlleles))
  else if(length(alleles) != length(oldAlleles))
    stop2("Replacement `alleles` has incorrect length (should be either 1 or ",length(oldAlleles), ")")

  if(!is.matrix(alleles))
    dim(alleles) = dim(oldAlleles)

  if(is.null(rownames(alleles)))
    rownames(alleles) = ids
  else if(!setequal(rownames(alleles), ids))
    stop2("Unknown ID label(s) found in rownames of `alleles`: ",
          setdiff(ids, rownames(alleles)))

  # Function for doing one component
  setAllelesComponent = function(comp) {
    ids_comp = intersect(ids, labels(comp))

    am = completeAlleleMatrix[labels(comp), , drop = FALSE]
    am[ids_comp, ] = alleles[ids_comp, ]

    # Locus attributes
    markers = markers %||% seq_markers(comp)
    loci = getLocusAttributes(comp, markers = markers)

    # Convert back to marker list and replace the old ones
    mlistNew = allelematrix2markerlist(comp, am, locusAttributes = loci, missing = NA)

    midx = whichMarkers(comp, markers)
    comp$MARKERS[midx] = mlistNew

    # Return modified ped
    comp
  }


  if(is.pedList(x))
    lapply(x, setAllelesComponent)
  else
    setAllelesComponent(x)
}


# For internal use
allelematrix2markerlist = function(x, alleleMatrix, locusAttributes, missing = 0, sep = NULL, validate = TRUE) {

  if(!is.matrix(alleleMatrix) && !is.data.frame(alleleMatrix))
    stop2("Argument `alleleMatrix` must be either a matrix or a data.frame")

  m = as.matrix(alleleMatrix)

  # Marker names in matrix (if present)
  hasMatrixNames = !is.null(nms <- colnames(m)) && !any(is.na(nms))

  # Rownames (ID labels)
  row_nms = rownames(m)

  # If rownames, reorder according to pedigree
  if(!is.null(row_nms)) {

    tmp = matrix("0", nrow = pedsize(x), ncol = ncol(m))
    idx = match(row_nms, labels(x))

    tmp[idx[!is.na(idx)], ] = m[row_nms[!is.na(idx)], ]

    m = tmp
  }
  else {  # If no rownames - dimensions must be correct
    if(nrow(m) != pedsize(x))
    stop2(sprintf("Incompatible input.\n  Pedigree size = %d\n  Allele matrix rows = %d",
                  pedsize(x), nrow(m)))
  }

  # If `sep` is given, split each column into single-allele columns
  if(!is.null(sep)) {
    m = split_genotype_cols(m, sep, missing)
  }
  else {
    if (ncol(m) %% 2 != 0)
      stop2("Uneven number of marker allele columns")

    # Marker names: Use odd numbered columns
    if(hasMatrixNames) {
      newnms = nms[seq.int(1, length(nms), by = 2)]

      # M1.1, M1.2, M2.1, M2.2, ... --> M1, M2, ...
      if(all(endsWith(nms, ".1") | endsWith(nms, ".2")))
        newnms = substr(newnms, 1, nchar(newnms) - 2)

      nms = newnms

      # Convert numerical names to NA NB: Skipping this caused problems in pedbuildr
      if (isTRUE(any(suppressWarnings(nms == as.integer(nms)))))
        hasMatrixNames = FALSE
    }
  }

  # Settle the number of markers
  nMark = ncol(m)/2

  # Check for (nontrivial) duplicated marker names found in matrix
  if(hasMatrixNames) {
    dups = duplicated(nms) & !is.na(nms) & nms != ""
    if(any(dups))
      stop2("Duplicated marker name: ", nms[dups])
  }

  # Replace `missing` with zeroes
  if(!identical(missing, 0))
    m[m %in% missing] = 0

  # Quick return if no locus attributes are given
  if(is.null(locusAttributes)) {
    mlist = lapply(seq_len(nMark), function(i) {
      mi = m[, c(2*i - 1, 2*i), drop = FALSE]
      nm = nms[i] # NULL is ok!
      marker(x, allelematrix = mi, name = nm, validate = FALSE)
    })

    return(mlist)
  }

  ####### locusAttributes #########

  # If same attributes for all: Recycle
  if (length(locusAttributes) == 1 && nMark > 1) {
    nm = locusAttributes[[1]]$name
    if(!is.null(nm) && !is.na(nm))
      stop2("Cannot recycle `name` attribute")
    locusAttributes = rep(locusAttributes, nMark)
  }

  # Scenario 1: Allele matrix has marker names
  if(hasMatrixNames) {

    # Use names(locusAttributes) if these exist
    nms_attr = names(locusAttributes)
    hasAttrNames = !is.null(nms_attr)

    # Otherwise use `name` attributes if given
    if(!hasAttrNames) {
      nms_attr = vapply(locusAttributes, function(a) as.character(a$name %||% NA), "")
      hasAttrNames = !all(is.na(nms_attr))
    }

    if(hasAttrNames) {
      if(anyNA(idx <- match(nms, nms_attr)))
        stop2("Marker name found in `allelematrix`, but not in `locusAttributes`: ", setdiff(nms, nms_attr))
      locusAttributes = locusAttributes[idx]
    }
    else {
      # If no names found in locusAttributes: Insert names from matrix!
      if(nMark != length(locusAttributes))
        stop2("When `locusAttributes` doesn't contain marker names, its length must match the number of markers")
      locusAttributes = lapply(seq_len(nMark), function(i) {
        a = locusAttributes[[i]]
        a$name = nms[i]
        a
      })
    }
  }
  else {
    # No marker names in matrix
    if(nMark != length(locusAttributes))
      stop2("Length of `locusAttributes` must match `alleleMatrix`, when the latter doesn't contain marker names")
  }

  mlist = lapply(seq_len(nMark), function(i) {
    attri = locusAttributes[[i]]
    attri$x = x
    attri$allelematrix = m[, c(2*i - 1, 2*i), drop = FALSE]
    attri$validateMut = FALSE
    do.call(marker, attri) # create marker object
  })

  class(mlist) = "markerList"
  mlist
}


split_genotype_cols = function(m, sep, missing) {
  nas = is.na(m) | m == missing
  if(all(nas))
    return(matrix(0, nrow = nrow(m), ncol = 2*ncol(m)))

  nonNA = m[!nas][1]
  if(!grepl(sep, nonNA))
    stop2("Allele separator not found in first non-NA entry of allele matrix: ", nonNA)

  # Replace NA's and missing by <miss>/<miss>. (Suboptimal strategy, but simple)
  m[nas] = sprintf("%s%s%s", missing, sep, missing)

  nc = ncol(m)
  nr = nrow(m)
  splitvec = unlist(strsplit(m, sep, fixed = TRUE))
  msplit = matrix(0, ncol = 2 * nc, nrow = nr)
  msplit[, 2 * seq_len(nc) - 1] = splitvec[2 * seq_len(nc * nr) - 1]
  msplit[, 2 * seq_len(nc)] = splitvec[2 * seq_len(nc * nr)]
  msplit
}

# For internal use
markerlist2allelematrix = function(mlist, missing = NA) {

  # List of vectors, each of length 2*pedsize
  allelelist = lapply(mlist, function(m) c(missing, alleles.marker(m))[m + 1])

  # Transform to matrix (faster than cbind)
  amat = unlist(allelelist)
  nMark = length(mlist)
  nInd = length(amat)/(2*nMark)
  dim(amat) = c(nInd, 2*nMark)

  # Column names
  mnames = sapply(mlist, name)
  if(any(naname <- is.na(mnames)))
    mnames[naname] = as.character(which(naname))
  colnames(amat) = paste0(rep(mnames, each = 2), c(".1", ".2"))

  # Return character matrix
  amat
}

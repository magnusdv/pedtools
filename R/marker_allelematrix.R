#' Allele matrix manipulation
#'
#' Functions for getting and setting the genotypes of multiple
#' individuals/markers simultaneously
#'
#' If the `alleles` argument of `setAlleles()` is not a matrix, it is recycled
#' (if neccessary), and converted into a matrix of the correct dimensions. For
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
#'   modified alleles. In particular, all locus annotations are unchanged.
#'
#' @seealso [transferMarkers()]
#'
#' @examples
#' x = nuclearPed(1)
#' m1 = marker(x, `2` = 1:2, alleles = 1:2, name="m1")
#' m2 = marker(x, `3` = 2, alleles = 1:2, name="m2")
#' x = setMarkers(x, list(m1, m2))
#'
#' mat1 = getAlleles(x)
#' mat2 = getAlleles(x, ids=2:3, markers = "m2")
#' stopifnot(identical(mat1[2:3, 3:4], mat2))
#'
#' # Remove all genotypes
#' y = setAlleles(x, alleles = 0)
#' y
#'
#' # Setting a single genotype
#' z = setAlleles(y, ids="1", marker = "m2", alleles = 1:2)
#'
#' # Alternative: In-place modification with `genotype()`
#' genotype(y, id = "1", marker = "m2") = 1:2
#' stopifnot(identical(y,z))
#'
#' ### Manipulation of pedlist objects
#' s = transferMarkers(x, singleton("s"))
#' peds = list(x, s)
#'
#' getAlleles(peds)
#'
#' setAlleles(peds, ids="s", marker="m1", alleles=1:2)
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
    if(!is.null(ids) && !all(ids %in% unlist(lapply(x, labels))))
      stop2("Unknown ID label: ", setdiff(ids, unlist(lapply(x, labels))))

    # Check equality of marker counts and names
    mNames = lapply(x, function(comp) name(comp, markers = seq_along(nMarkers(comp))))
    if(length(unique(mNames)) > 1)
      stop2("Components cannot have different number of markers, or different marker names. Please file an issue if this is important to you.")

    # Extract alleles from each component
    amList = lapply(x, function(comp) {
      ids_comp = if(!is.null(ids)) intersect(ids, labels(comp))
      getAlleles(comp, ids = ids_comp, markers = markers)
    })

    # Bind into single matrix
    res = do.call(rbind, amList)

    # Sort rows according to input `ids`
    if(!is.null(ids))
      res = res[as.character(ids), , drop = F]

    return(res)
  }

  # Main body: x is now is single `ped` object
  if(is.null(ids))
    ids = labels(x)

  if(length(ids) == 0)
    return(NULL)

  if(is.null(markers))
    markers = seq_len(nMarkers(x))

  mlist = getMarkers(x, markers = markers)
  if(length(mlist) == 0)
    return(NULL)

  # Convert to character matrix
  am = markerlist2allelematrix(mlist)

  # Subset and set rownames
  amSubset = am[internalID(x, ids), , drop = F]
  rownames(amSubset) = ids

  # Return
  amSubset
}

#' @rdname getAlleles
#' @export
setAlleles = function(x, ids = NULL, markers = NULL, alleles) {
  if(!is.ped(x) && !is.pedList(x))
    stop2("The first argument must be a `ped` object or a list of such")

  if(is.pedList(x)) {

    if(is.matrix(alleles) && is.null(rownames(alleles))) {
      if(is.null(ids))
        rownames(alleles) = unlist(lapply(x, labels))
      else
        stop2("When `ids` is non-NULL and `alleles` is a matrix, it must have rownames")
    }

    # Set alleles one component at a time
    y = lapply(x, function(comp) {
      ids_comp = if(is.null(ids)) labels(comp) else intersect(ids, labels(comp))
      am_comp = if(is.matrix(alleles)) alleles[ids_comp, ] else alleles
      setAlleles(comp, ids = ids_comp, markers = markers, alleles = am_comp)
    })
    return(y)
  }

  # Main body: x is now is single `ped` object
  if(is.null(ids))
    ids = labels(x)

  if(length(ids) == 0)
    return(x)

  if(is.null(markers))
    markers = seq_len(nMarkers(x))

  midx = whichMarkers(x, markers = markers)
  if(length(midx) == 0) {
    cat("No markers to set, returning `x` unchanged\n")
    return(x)
  }

  # Fix `alleles` if not matrix
  if(is.data.frame(alleles))
    alleles = as.matrix(alleles)
  if(!is.matrix(alleles))
    alleles = matrix(alleles, nrow = length(ids), ncol = 2*length(markers))

  mlist = getMarkers(x, midx)
  am = markerlist2allelematrix(mlist)
  am[internalID(x, ids), ] = alleles

  loci = lapply(mlist, attributes)
  mlistNew = allelematrix2markerlist(x, am, locus_annotations = loci, missing = NA)

  x$markerdata[midx] = mlistNew
  x
}


# For internal use
allelematrix2markerlist = function(x, allele_matrix, locus_annotations, missing=0, allele_sep=NULL) {

  if(!is.matrix(allele_matrix) && !is.data.frame(allele_matrix))
    stop2("Argument `allele_matrix` must be either a matrix or a data.frame")

  m = as.matrix(allele_matrix)
  row_nms = rownames(m)

  # If no rownames - dimensions must be correct
  if(is.null(row_nms)) {
    if(nrow(m) != pedsize(x))
    stop2("Incompatible input.\n  Pedigree size = ", pedsize(x),
          "\n  Allele matrix rows = ", nrow(m))
  }
  else {
    tmp = matrix("0", nrow = pedsize(x), ncol = ncol(m))
    idx = match(row_nms, labels(x))

    tmp[idx[!is.na(idx)], ] = m[row_nms[!is.na(idx)], ]

    m = tmp
    #   missing_labs = setdiff(labels(x), row_nms)
    #   if (length(missing_labs))
    #     stop2("Pedigree member missing in allele matrix row names: ", missing_labs)
    #   unknown_labs = setdiff(row_nms, labels(x))
    #   if (length(unknown_labs))
    #     stop2("Unknown row name in allele matrix: ", unknown_labs)
    #
    #   # Reorder
    #   if(!identical(labels(x), row_nms))
    #     m = m[labels(x), ]
    #
  }

  # If allele_sep is given, interpret entries as diploid genotypes
  if(!is.null(allele_sep))
    m = split_genotype_cols(m, allele_sep, missing)

  if (ncol(m) %% 2 != 0)
    stop2("Uneven number of marker allele columns")

  if(!identical(missing, 0)) {
    m[m %in% missing] = 0
  }

  nMark = ncol(m)/2
  ann = locus_annotations

  # Quick return if no annotations given
  if(is.null(ann)) {
    mlist = lapply(seq_len(nMark), function(i) {
      mi = m[, c(2*i - 1, 2*i), drop = FALSE]
      marker(x, allelematrix = mi, validate = FALSE)
    })

    return(mlist)
  }

  # If same annotations for all: Recycle
  if (!is.list(ann[[1]]))
    ann = rep(list(ann), nMark)

  if (length(ann) != nMark)
    stop2(sprintf("Length of annotation list (%d) does not equal number of markers (%d)",
                  length(ann), nMark))

  mlist = lapply(seq_len(nMark), function(i) {
    attribs = ann[[i]]
    attribs$x = x
    attribs$allelematrix = m[, c(2*i - 1, 2*i), drop = FALSE]
    do.call(marker, attribs)
  })

  class(mlist) = "markerList"
  mlist
}

split_genotype_cols = function(m, allele_sep, missing) {
  nas = is.na(m) | m == missing
  if(all(nas))
    return(matrix(0, nrow = nrow(m), ncol = 2*ncol(m)))

  nonNA = m[!nas][1]
  if(!grepl(allele_sep, nonNA))
    stop2("Allele separator not found in first non-NA entry of allele matrix: ", nonNA)

  # Replace NA's and missing by <miss>/<miss>. (Suboptimal strategy, but simple)
  m[nas] = sprintf("%s%s%s", missing, allele_sep, missing)

  nc = ncol(m)
  nr = nrow(m)
  splitvec = unlist(strsplit(m, allele_sep, fixed = T))
  msplit = matrix(0, ncol = 2 * nc, nrow = nr)
  msplit[, 2 * seq_len(nc) - 1] = splitvec[2 * seq_len(nc * nr) - 1]
  msplit[, 2 * seq_len(nc)] = splitvec[2 * seq_len(nc * nr)]
  msplit
}

# For internal use
markerlist2allelematrix = function(mlist, missing = NA) {
  # List of 2-col character matrices
  alist = lapply(mlist, function(m) {
    a = c(missing, alleles(m))[m + 1]
    dim(a) = dim(m)
    a
  })

  # Bind
  amat = do.call(cbind, alist)

  # Column names
  mnames = sapply(mlist, name)
  if(any(naname <- is.na(mnames)))
    mnames[naname] = paste0("na", 1:sum(naname))
  colnames(amat) = paste(rep(mnames, each=2), 1:2, sep=".")

  # Return character matrix
  amat
}

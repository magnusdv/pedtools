#' Attach markers to pedigrees
#'
#' In many applications it is useful to _attach_ markers to their associated
#' `ped` object. In particular for bigger projects with many markers, this makes
#' it easier to manipulate the dataset as a unit. The function `setMarkers()`
#' replaces all existing markers with the supplied ones, while `addMarkers()`
#' appends the supplied markers to any existing ones.
#'
#' @param x A `ped` object
#' @param m Either a single `marker` object or a list of `marker` objects
#' @param allele_matrix A matrix with `pedsize(x)` rows, containing the observed
#'   alleles for one or several markers. The matrix must have either 1 or 2
#'   columns per marker. If the former, then the argument `allele_sep` must be
#'   non-NULL, and will be used to split all entries.
#' @param locus_annotations A list of lists, with annotations for each marker.
#'   See [marker()] for possible entries.
#' @param missing A single character (or coercible to one) indicating the symbol
#'   for missing alleles.
#' @param allele_sep If this is a single string, each entry of `allele_matrix`
#'   is interpreted as a genotype, and will be split by calling `str_split(...,
#'   split = allele_sep, fixed = T)`. For example, if the entries are formatted
#'   as "A/B", put `allele_sep="/"`. Default: NULL.
#'
#' @return A `ped` object.
#' @examples
#' x = singleton(1)
#' m1 = marker(x, '1'=1:2)
#' m2 = marker(x, '1'='a')
#'
#' x = setMarkers(x, m1)
#' x = addMarkers(x, m2)
#' x
#'
#' # Reversing the order of the markers
#' x = setMarkers(x, list(m2, m1))
#' x
#'
#' @name marker_attach
NULL

#' @rdname marker_attach
#' @export
setMarkers = function(x, m = NULL, allele_matrix = NULL, locus_annotations = NULL, missing=0, allele_sep=NULL) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(is.null(m) && is.null(allele_matrix)) {
    x['markerdata'] = list(NULL)
    return(x)
  }

  if (is.marker(m))
    mlist = list(m)
  else if (is.markerList(m))
    mlist = m
  else if (is.null(m))
    mlist = allelematrix2markerlist(x, allele_matrix, locus_annotations, missing, allele_sep)
  else
    stop2("Argument `m` must be either a single `marker` object, a list of such, or NULL")

  class(mlist) = "markerList"
  checkConsistency(x, mlist)
  x$markerdata = mlist
  x
}

#' @rdname marker_attach
#' @export
addMarkers = function(x, m = NULL, allele_matrix = NULL, locus_annotations = NULL, missing=0, allele_sep=NULL) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  if (is.marker(m))
    mlist = list(m)
  else if (is.markerList(m))
    mlist = m
  else if (is.null(m))
    mlist = allelematrix2markerlist(x, allele_matrix, locus_annotations, missing, allele_sep)
  else
    stop2("Argument `m` must be either a single `marker` object, a list of such, or NULL")

  # Append to x and return
  checkConsistency(x, mlist)
  mlist = c(x$markerdata, mlist)
  class(mlist) = "markerList"
  x$markerdata = mlist
  x
}



allelematrix2markerlist = function(x, allele_matrix, locus_annotations, missing=0, allele_sep=NULL) {

  if(!is.matrix(allele_matrix) && !is.data.frame(allele_matrix))
    stop2("Argument `allele_matrix` must be either a matrix or a data.frame")

  m = as.matrix(allele_matrix)
  if(nrow(m) != pedsize(x))
    stop2("Incompatible input.\n  Pedigree size = ", pedsize(x),
         "\n  Allele matrix rows = ", nrow(m))

  # # If row names are given, use them to re-order matrix
  # if (!is.null(row_nms <- rownames(m))) {
  #
  #   # Check compatibility
  #   missing_labs = setdiff(x$LABELS, row_nms)
  #   if (length(missing_labs))
  #     stop2("Pedigree member missing in allele matrix row names: ", missing_labs)
  #   unknown_labs = setdiff(row_nms, x$LABELS)
  #   if (length(unknown_labs))
  #     stop2("Unknown row name in allele matrix: ", unknown_labs)
  #
  #   # Reorder
  #   if(!identical(x$LABELS, row_nms))
  #     m = m[x$LABELS, ]
  # }

  # If allele_sep is given, interpret each column as a marker
  if(!is.null(allele_sep))
    m = split_genotype_cols(m, allele_sep, missing)

  if (ncol(m) %% 2 != 0)
    stop2("Uneven number of marker allele columns")

  if(missing != 0) {
    m[m == missing] = 0
  }

  nMark = ncol(m)/2
  LABELS = x$LABELS
  SEX = x$SEX

  ann = locus_annotations

  # Quick return if no annotations given
  if(is.null(ann)) {
    mlist = lapply(seq_len(nMark), function(i) {
      mi = m[, c(2*i - 1, 2*i), drop = FALSE]
      .createMarkerObject(mi, pedmembers=LABELS, sex=SEX)
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
    if (is.null(attribs))
      return(NULL)

    attribs$matr = m[, c(2*i - 1, 2*i), drop = FALSE]
    attribs$pedmembers = LABELS
    attribs$sex = SEX
    do.call(.createMarkerObject, attribs)
  })

  mlist[unlist(lapply(mlist, is.null))] = NULL

  class(mlist) = "markerList"
  mlist
}

split_genotype_cols = function(m, allele_sep, missing) {
  nas = is.na(m)
  if(all(nas))
    return(matrix(0, nrow = nrow(m), ncol = 2*ncol(m)))

  nonNA = m[!nas][1]
  if(!grepl(allele_sep, nonNA))
    stop2("Allele separator not found in first non-NA entry of allele matrix: ", nonNA)

  # Replace NA's by <miss>/<miss>. (Suboptimal strategy, but simple)
  m[nas] = sprintf("%s%s%s", missing, allele_sep, missing)

  nc = ncol(m)
  nr = nrow(m)
  splitvec = unlist(strsplit(m, allele_sep, fixed = T))
  msplit = matrix(0, ncol = 2 * nc, nrow = nr)
  msplit[, 2 * seq_len(nc) - 1] = splitvec[2 * seq_len(nc * nr) - 1]
  msplit[, 2 * seq_len(nc)] = splitvec[2 * seq_len(nc * nr)]
  msplit
}

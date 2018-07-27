#' Attach marker objects to a pedigree
#'
#' @param x A `ped` object
#' @param m Either a single `marker` object, a list of `marker` objects, or a data.frame or matrix.
#' @param annotations A list of marker annotations.
#'
#' @export
setMarkers = function(x, m, annotations = NULL) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  if (is.null(m)) {
    x['markerdata'] = list(NULL)
    return(x)
  }

  mlist = NULL
  if (is.marker(m))
    mlist = list(m)
  else if (is.markerList(m))
    mlist = m
  else if (!(is.data.frame(m) || is.matrix(m)))
    stop("Argument 'm' must be either:\n",
         " * a single `marker` object`\n",
         " * a list of `marker` objects\n",
         " * a data.frame or matrix.", call.=F)

  # If markerlist, attach to x and return
  if(!is.null(mlist)) {
    class(mlist) = "markerList"
    checkConsistency(x, mlist)
    x$markerdata = mlist
    return(x)
  }

  m = as.matrix(m)
  nc = ncol(m)
  nr = nrow(m)
  if(nr != pedsize(x))
    stop("Incompatible input. Pedigree has size ", pedsize(x),
         " but allele matrix has ", nr, " rows.", calls.=FALSE)

  # If input has 1 genotype per column, split alleles to separate columns
  is_merged = is.character(m) && grepl("/", m[1,1])
  if (is_merged) {
    splitvec = unlist(strsplit(m, "/", fixed = T))
    msplit = matrix(0, ncol = 2 * nc, nrow = nr)
    msplit[, 2 * seq_len(nc) - 1] = splitvec[2 * seq_len(nc * nr) - 1]
    msplit[, 2 * seq_len(nc)] = splitvec[2 * seq_len(nc * nr)]
    m = msplit
  }

  if (ncol(m) %% 2 != 0)
    stop2("Uneven number of marker allele columns")

  nMark = ncol(m)/2

  if (!is.null(annotations)) {
    if (length(annotations) == 2 && !is.null(names(annotations)))
      annotations = rep(list(annotations), nMark)  # if given attrs for a single marker
    else if (length(annotations) != nMark)
      stop2(sprintf("Length of annotation list (%d) does not equal number of markers (%d)",
                    length(annotations), nMark))

    mlist = lapply(1:nMark, function(i) {
      if (is.null(attribs <- annotations[[i]]))
        return(NULL)
      mi = m[, c(2 * i - 1, 2 * i), drop = FALSE]
      do.call(.createMarkerObject, c(list(matr = mi), attribs))
    })
  }
  else {
    mlist = lapply(1:nMark, function(i) {
      mi = m[, c(2 * i - 1, 2 * i), drop = FALSE]
      .createMarkerObject(mi)
    })
  }
  mlist[sapply(mlist, is.null)] = NULL
  class(mlist) = "markerList"
  x$markerdata = mlist
  x
}

#' @rdname setMarkers
#' @export
addMarkers = function(x, m, annotations = NULL) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  if (is.null(m))
    return(x)

  mlist = NULL
  if (is.marker(m))
    mlist = list(m)
  else if (is.markerList(m))
    mlist = m
  else stop2("Matrix or data.frame input is not implemented yet")

  # If markerlist, append to x and return
  if(!is.null(mlist)) {
    checkConsistency(x, mlist)
    x$markerdata = c(x$markerdata, mlist)
    class(x$markerdata) = "markerList"
    return(x)
  }
}

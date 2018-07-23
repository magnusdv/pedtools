#' Attach marker objects to a pedigree
#'
#' @param x A `ped` object
#' @param m Either a single `marker` object, a list of `marker` objects, or a data.frame or matrix.
#' @param annotations A list of marker annotations.
#'
#' @export
setMarkers = function(x, m, annotations = NULL) {
  assert_that(is.ped(x))

  if (is.null(m))
    mlist = NULL
  else if (is.marker(m))
    mlist = list(m)
  else if (is.markerList(m))
    mlist = m
  else if (!(is.data.frame(m) || is.matrix(m)))
    stop("Argument 'm' must be either:\n",
         " * a single `marker` object`\n",
         " * a list of `marker` objects\n",
         " * a data.frame or matrix.", call.=F)

  # If any of the above kicked in, append to x and return
  if(exists("mlist", environment(), inherits=F)) {
    class(mlist) = "markerList"
    checkConsistency(x, mlist)
    x$markerdata = mlist
    return(x)
  }

  m = as.matrix(m)
  nc = ncol(m)
  nr = nrow(m)
  if(nr != pedSize(x))
    stop("Incompatible input. Pedigree has size ", pedSize(x),
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
    stop("Uneven number of marker allele columns")

  nMark = ncol(m)/2

  if (!is.null(annotations)) {
    if (length(annotations) == 2 && !is.null(names(annotations)))
      annotations = rep(list(annotations), nMark)  # if given attrs for a single marker
    else if (length(annotations) != nMark)
      stop("Length of marker annotation list does not equal number of markers.")

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

#' @export
addMarker = function(x, m, ...) {
  N = pedSize(x)
  if (is.matrix(m) || is.data.frame(m))
    stopifnot(nrow(m) == N, ncol(m) == 2)
  if (inherits(m, "marker"))
    m = list(m)
  if (!is.list(m) && length(m) == 1)
    m = matrix(m, ncol = 2, nrow = N)  #gives a nice way of setting an empty or everywhere-homozygous marker, e.g.: x=addMarker(x,0)
  mm = .createMarkerObject(m, ...)
  setMarkers(x, structure(c(x$markerdata, list(mm)), class = "markerdata"))
}

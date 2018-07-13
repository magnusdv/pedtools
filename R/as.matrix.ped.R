#' Convert ped to matrix
#'
#' Converts a `ped` object to a numeric matrix using internal labels,
#' with additional info neccessary to recreate the original `ped`
#' attached as attributes.
#'
#' `restore_ped` is the reverse of `as.matrix.ped`.
#'
#' @param x a [ped()] object. In `restore_ped`: A
#' numerical matrix.
#' @param include.attrs a logical indicating if marker annotations and other
#' info should be attached as attributes. See value.
#' @param attrs a list containing labels and other `ped`
#' info compatible with `x`, in the format produced by `as.matrix`.
#' If NULL, the attributes of `x` itself are used.
#' @param check a logical, forwarded to [ped()]. If FALSE, no
#' checks for pedigree errors are performed.
#' @param \dots not used.
#'
#' @return For `as.matrix`: A matrix with `x$NIND` rows.
#' If `include.attrs = TRUE` the matrix has the following attributes:
#' \itemize{
#' \item{`labels`}{ a list of marker annotations}
#' \item{`famid`}{ the availability vector}
#' }
#'
#' For `restore_ped`: A `ped` object.
#' @author Magnus Dehli Vigeland
#' @seealso [ped()]
#'
#' @examples
#'
#' x = setLabels(nuclearPed(1), letters[1:3])
#'
#' # To examplify the ped -> matrix -> ped trick, we show how to
#' # reverse the internal ordering of the pedigree.
#' m = as.matrix(x, include.attrs=TRUE)
#' m[] = m[3:1, ]
#'
#' # Must reverse the labels also:
#' attrs = attributes(m)
#' attrs$labels = rev(attrs$labels)
#'
#' # Restore ped:
#' y = restore_ped(m, attrs=attrs)
#'
#' # Of course a simpler way is use reorder():
#' z = reorder(x, 3:1)
#' stopifnot(identical(y, z))
#'
#' @export
as.matrix.ped = function(x, include.attrs = TRUE, ...) {
  m = cbind(x$ID, x$FID, x$MID, x$SEX)
  if(hasMarkers(x)) {
    markermatr = do.call(cbind, x$markerdata)
    m = cbind(m, markermatr)
  }
  if (include.attrs) {
    attr(m, "famid") = x$FAMID
    attr(m, "labels") = x$LABELS
    attr(m, "hasLoops") = x$hasLoops
    attr(m, "loop_breakers") = x$loop_breakers
    if(hasMarkers(x)) {
      attr(m, "markerattr") = lapply(x$markerdata, attributes)
    }
  }
  m
}


#' @rdname as.matrix.ped
#' @export
restore_ped = function(x, attrs = NULL, check = TRUE) {
  if (is.null(attrs))
    attrs = attributes(x)
  p = ped(id=x[,1], fid=x[,2], mid=x[,3], sex=x[,4], famid=attr(x, "famid"), check = check, reorder=F)
  p = setLabels(p, attrs$labels)
  p['loop_breakers'] = list(attrs$loop_breakers) # Trick to keep explicit NULLs

  ### Markers
  if((nc <- ncol(x)) > 4) {
    if(nc %% 2 != 0) stop("The number of allele columns must be even")
    markerattr = attrs$markerattr

    mlist = lapply(seq_len((nc-4)/2), function(k) {
      m = x[, c(3 + 2*k, 4 + 2*k), drop = F]
      attr = markerattr[[k]]
      attr$dim = dim(m)
      attributes(m) = attr
      m
    })
    class(mlist) = "markerList"
    p = setMarkers(p, mlist)
  }
  p
}

#' Convert ped to data.frame
#'
#' Convert a `ped` object to a data.frame with columns id, fid, mid and sex. The
#' output uses the original ID labels, unlike [as.matrix.ped()] which uses
#' internal numeric IDs. The latter is safer for manipulating the pedigree
#' structure, since the internal IDs are always distinct. (See [setLabels()].)
#' The main use of the data.frame method is printing.
#'
#' @param x Object of class `ped`.
#' @param ... Further parameters
#' @param markers (Optional) Vector of marker indices. By default, all markers
#'   are included.
#'
#' @export
as.data.frame.ped = function(x, ..., markers) {
  lab = x$LABELS
  fid = mid = rep("0", pedSize(x))
  fid[x$FID > 0] = lab[x$FID]
  mid[x$MID > 0] = lab[x$MID]
  df = data.frame(id = lab, fid=fid, mid=mid, sex=x$SEX, stringsAsFactors = F)

  if(hasMarkers(x)) {
    mlist = if(missing(markers)) x$markerdata else getMarkers(x, markers)
    geno = .prettyMarkers(mlist, ...)
    df = cbind(df, geno)
  }
  df
}

#' Printing pedigrees
#'
#' Print a `ped` object using original labels.
#'
#' This first calls [as.data.frame.ped()] and then prints the resulting
#' data.frame. The data.frame is returned invisibly.
#'
#' @param x object of class `ped`.
#' @param ... (optional) arguments passed on to [print.data.frame()].
#' @param markers (optional) vector of marker indices. If missing, and `x` has
#'   less than 10 markers, they are all displayed. If `x` has 10 or more
#'   markers, the first 5 are displayed.
#' @param verbose If TRUE, a message is printed if only the first 5 markers are
#'   printed. (See above).
#' @export
print.ped = function(x, ..., markers, verbose=TRUE) {
  nm = nMarkers(x)
  showmess = F
  if (missing(markers)) {
    if (nm < 10)
      markers = seq_len(nm)
    else {
      markers = 1:5
      showmess = T
    }
  }
  else {
    if (any(markers > nm)) stop("Nonexisting marker(s) indicated")
  }
  datafr = as.data.frame(x, markers=markers)
  print(datafr, ...)

  if(showmess && verbose)
    message("Only 5 (out of ", nm, ") markers are shown. See `?print.ped` for options.")

  invisible(datafr)
}


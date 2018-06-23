#' Convert ped to matrix
#'
#' Converts a \code{ped} object to a numeric matrix using internal labels,
#' with additional info neccessary to recreate the original \code{ped}
#' attached as attributes.
#'
#' \code{restore_ped} is the reverse of \code{as.matrix.ped}.
#'
#' @param x a \code{\link{ped}} object. In \code{restore_ped}: A
#' numerical matrix.
#' @param include.attrs a logical indicating if marker annotations and other
#' info should be attached as attributes. See value.
#' @param attrs a list containing labels and other \code{ped}
#' info compatible with \code{x}, in the format produced by \code{as.matrix}.
#' If NULL, the attributes of \code{x} itself are used.
#' @param check a logical, forwarded to \code{\link{ped}}. If FALSE, no
#' checks for pedigree errors are performed.
#' @param \dots not used.
#'
#' @return For \code{as.matrix}: A matrix with \code{x$NIND} rows.
#' If \code{include.attrs = TRUE} the matrix has the following attributes:
#' \itemize{
#' \item{\code{labels}}{ a list of marker annotations}
#' \item{\code{famid}}{ the availability vector}
#' }
#'
#' For \code{restore_ped}: A \code{ped} object.
#' @author Magnus Dehli Vigeland
#' @seealso \code{\link{ped}}
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
    if (include.attrs) {
        attr(m, "famid") = x$FAMID
        attr(m, "labels") = x$LABELS
        attr(m, "hasLoops") = x$hasLoops
        attr(m, "loop_breakers") = x$loop_breakers
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
  p
}

#' Convert ped to data.frame
#'
#' Convert a \code{ped} object to a data.frame with columns id, fid, mid and sex.
#' The output uses the original ID labels, unlike \code{\link{as.matrix.ped}} which uses
#' internal numeric IDs. The latter is safer for manipulating the pedigree structure,
#' since the internal IDs are always distinct. (See \code{\link{setLabels}}.)
#' The main use of the data.frame method is printing.
#'
#' @param x object of class \code{ped}.
#' @param ... Not used.
#'
#' @export
as.data.frame.ped = function(x, ...) {
  lab = x$LABELS
  fid = mid = rep("0", pedSize(x))
  fid[x$FID > 0] = lab[x$FID]
  mid[x$MID > 0] = lab[x$MID]
  data.frame(id = lab, fid=fid, mid=mid, sex=x$SEX, stringsAsFactors = F)
}

#' Printing pedigrees
#'
#' Print a \code{ped} object using original labels.
#'
#' This first calls \code{\link{as.data.frame.ped}} and then prints the
#' resulting data.frame. The data.frame is returned invisibly.
#'
#' @param x object of class \code{ped}.
#' @param ... optional arguments passed on to \code{\link{print.data.frame}}.
#'
#' @export
print.ped = function(x, ...) {
  datafr = as.data.frame(x)
  print(datafr, ...)
  invisible(datafr)
}


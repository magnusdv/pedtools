#' Merge two pedigrees
#'
#' This function merges two ped objects, joining them at the individuals
#' with equal ID labels. This is especially useful for building 'top-heavy'
#' pedigrees. Only ped objects without marker data are supported.
#'
#'
#' @param x,y \code{\link{ped}} objects
#' @param ... further arguments passed along to \code{\link{ped}}, e.g.
#'   \code{famid}, \code{check} and \code{reorder}.
#'
#' @return A \code{ped} object.
#' @author Magnus Dehli Vigeland
#'
#' @examples
#'
#' # Creating a trio where each parent have first cousin parents.
#' # (Alternatively, this could be built using many calls to addParents().)
#'
#' x = cousinsPed(1)
#' x = addChildren(x, father=5, mother=8, nch=1, ids=9)
#' x = addChildren(x, father=9, mother=10, nch=1, ids=11)
#'
#' y = relabel(cousinsPed(1), 101:108)
#' y = addChildren(y, father=105, mother=108, nch=1, sex=2, id=10)
#' y = addChildren(y, father=9, mother=10, nch=1, id=11)
#'
#' # Joining x and y at the common individuals 9,10,11:
#' z = mergePed(x,y)
#'
#' # plot all three pedigrees
#' par(mfrow=c(1,3))
#' plot(x); plot(y); plot(z)
#'
#' @export
mergePed = function(x, y, ...) {
  if (!is.null(x$markerdata) || !is.null(y$markerdata))
      stop("Merging is only supported for pedigrees without marker data", call.=FALSE)
  ids = intersect(x$LABELS, y$LABELS)
  if (length(ids) == 0)
      stop("Merging impossible: No common IDs", call.=FALSE)

  del = list(x = numeric(), y = numeric())
  for (i in ids) {
      if (getSex(x, i) != getSex(y, i))
          stop("Gender mismatch for individual ", i, call.=FALSE)
      parx = parents(x, i)
      pary = parents(y, i)
      if (length(pary) == 0)
          del$y = c(del$y, i) else if (length(parx) == 0)
          del$x = c(del$x, i) else if (all(parx == pary))
          del$y = c(del$y, i) else stop("Parent mismatch for individual ", i, call.=FALSE)
  }
  xm = as.data.frame(x)
  ym = as.data.frame(y)
  z = rbind(xm[!x$LABELS %in% del$x, ], ym[!y$LABELS %in% del$y, ])

  ped(id=z$id, fid=z$fid, mid=z$mid, sex=z$sex, ...)
}

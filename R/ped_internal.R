#' Internal ordering of pedigree members
#'
#' These functions give access to - and enable modifications of - the order in
#' which the members of a pedigree are stored. (This is the order in which the
#' members are listed when a `ped` object is printed to the screen.)
#'
#' The internal ordering is usually of little importance for end users, with one
#' important exception: Certain pedigree-traversing algorithms require parents
#' to precede their children. A special function, `parents_before_children()` is
#' provided for this purpose. This is a wrapper of the more general
#' `reorderPed()` which allows any permutation of the members.
#'
#' It should be noted that [ped()] by default calls `parents_before_children()`
#' whenever a pedigree is created, unless explicitly avoided with
#' `reorder=FALSE`.
#'
#' `has_parents_before_children()` can be used as a quick test to decide if it
#' is neccessary to call `parents_before_children()`.
#'
#' The utility `internalID()` converts ID labels to indices in the internal
#' ordering.
#'
#' @param x A `ped` object
#' @param neworder A permutation of `labels(x)` or of vector `1:pedsize(x)`. By
#'   default, the sorting order of the ID labels is used.
#' @param ids A character vector (or coercible to one) of original ID labels.
#'
#'
#' @seealso [ped()]
#' @examples
#' x = ped(id = 3:1, fid = c(1,0,0), mid = c(2,0,0), sex = c(1,2,1), reorder = FALSE)
#' x
#'
#' # The 'ids' argument is converted to character
#' internalID(x, ids = 3)
#' internalID(x, ids = "3")
#'
#' y = parents_before_children(x)
#' internalID(y, ids = 3)
#'
#' # A different ordering
#' reorderPed(x, c(2,1,3))
#'
#' @name ped_internal
NULL

#' @rdname ped_internal
#' @export
reorderPed = function(x, neworder = order(labels(x))) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(is.singleton(x))
    return(x)
  if(length(neworder) != pedsize(x))
    stop2("`neworder` must have length ", pedsize(x), ", not ", length(neworder))
  if(!setequal(neworder, labels(x)) && !setequal(neworder, 1:pedsize(x)))
    stop2("`neworder` must be a permutation of either `labels(x)` or `1:pedsize(x)`: ", neworder)
  if(is.character(neworder))
    neworder = internalID(x, neworder)

  xmatr = as.matrix(x)
  attr = attributes(xmatr)
  attr$LABELS = attr$LABELS[neworder]
  if(!is.null(lp <- attr$LOOP_BREAKERS))
    attr$LOOP_BREAKERS = matrix(match(lp, neworder), ncol=2)
  restore_ped(xmatr[neworder, ], attrs = attr)
}

#' @rdname ped_internal
#' @export
parents_before_children = function(x) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(is.singleton(x) || has_parents_before_children(x))
    return(x)

  neworder = 1:pedsize(x)
  i=1
  while (i < pedsize(x)) {
    current = neworder[i]
    maxpar = max(match(c(x$FIDX[current], x$MIDX[current]), neworder, nomatch = 0))
    if (maxpar > i) { # push current indiv below below parents
      neworder[i:maxpar] = neworder[c((i+1):maxpar, i)]
    }
    else i = i + 1
  }
  reorderPed(x, neworder)
}

#' @rdname ped_internal
#' @export
has_parents_before_children = function(x) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  idx = 1:pedsize(x)
  father_before_child = x$FIDX < idx
  mother_before_child = x$MIDX < idx
  all(father_before_child & mother_before_child)
}


#' @rdname ped_internal
#' @export
internalID = function(x, ids) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  int_ids = match(ids, labels(x))
  if (anyNA(int_ids))
    stop2("Unknown ID label: ", ids[is.na(int_ids)])
  int_ids
}

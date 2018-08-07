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
#' @param x a `ped` object
#' @param neworder a permutation of the vector `1:pedsize(x)`. By default, the
#'   order of the ID labels is used.
#' @param labels A character vector (or coercible to one) of original ID labels.
#'
#'
#' @seealso [ped()]
#' @examples
#' x = ped(3:1, fid=c(1,0,0), mid=c(2,0,0), sex=c(1,2,1), reorder=FALSE)
#' x
#'
#' # Note that 'label' is converted to character
#' internalID(x, label=3)
#'
#' y = parents_before_children(x)
#' internalID(y, label=3)
#'
#' # A different ordering
#' reorderPed(x, c(2,1,3))
#'
#' @name ped_internal
NULL

#' @rdname ped_internal
#' @export
reorderPed = function(x, neworder = order(x$LABELS)) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(is.singleton(x))
    return(x)
  xmatr = as.matrix(x)
  attr = attributes(xmatr)
  attr$labels = attr$labels[neworder]
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

  neworder = x$ID
  i=1
  while (i < pedsize(x)) {
    current = neworder[i]
    maxpar = max(match(c(x$FID[current], x$MID[current]), neworder, nomatch = 0))
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
  father_before_child = x$FID < x$ID
  mother_before_child = x$MID < x$ID
  all(father_before_child & mother_before_child)
}


#' @rdname ped_internal
#' @export
internalID = function(x, labels) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  int_ids = match(labels, x$LABELS)
  if (anyNA(int_ids))
    stop2("Unknown ID label: ", labels[is.na(int_ids)])
  int_ids
}

#' Internal ordering of pedigree members
#'
#' These functions give access to - and enable modifications of - the order in
#' which the members of a pedigree are stored. (This is the order in which the
#' members are listed when a `ped` object is printed to the screen.)
#'
#' While the internal pedigree ordering rarely matters, it is occasionally
#' important. The function `reorderPed()` permutes the internal ordering as
#' specified by the user. The most common use of this function is perhaps in
#' `parentsBeforeChildren()`, which ensures that all parents precede their
#' children. This is required by many pedigree-traversing algorithms.
#'
#' It should be noted that [ped()] by default calls `parentsBeforeChildren()`
#' whenever a pedigree is created, unless explicitly avoided with `reorder =
#' FALSE`.
#'
#' `hasParentsBeforeChildren()` can be used as a quick test to decide if it is
#' necessary to call `parentsBeforeChildren()`.
#'
#' The `foundersFirst()` function reorders the pedigree so that all the founders
#' come first.
#'
#' The utility `internalID()` converts ID labels to indices in the internal
#' ordering. If `x` is a list of pedigrees, the output is a data frame
#' containing both the component number and internal ID (within the component).
#'
#' @param x A `ped` object. Most of these functions also accepts ped lists.
#' @param neworder A permutation of `labels(x)` (or a subset of this),
#'   indicating the new internal ordering. If `internal = TRUE`, `neworder`
#'   refers to the internal ordering, so must be numeric. By default, the
#'   natural order of the ID labels is used.
#' @param internal A logical (default: FALSE). If TRUE, `neworder` is
#'   interpreted as referring to the internal ordering.
#' @param ids A character vector (or coercible to one) of original ID labels.
#' @param errorIfUnknown A logical. If TRUE (default), the function stops with
#'   an error if not all elements of `ids` are recognised as names of members in
#'   `x`.
#'
#' @seealso [ped()]
#'
#' @examples
#' x = ped(id = 3:1, fid = c(1,0,0), mid = c(2,0,0), sex = c(1,2,1), reorder = FALSE)
#' x
#'
#' # The 'ids' argument is converted to character, hence these are the same:
#' internalID(x, ids = 3)
#' internalID(x, ids = "3")
#'
#' hasParentsBeforeChildren(x)
#'
#' # Put parents first
#' parentsBeforeChildren(x)
#'
#' # Typical use of reorderPed: Swap sibling plot order
#' y = nuclearPed(2) |> reorderPed(4:3)
#' plot(y)
#'
#'
#' ### If labels are numeric, argument `internal` is important
#' z = singleton(1) |> addParents(1)
#' z
#' reorderPed(z, 1:3, internal = FALSE) # ID order = "1","2","3"
#' reorderPed(z, 1:3, internal = TRUE)  # index order: 1,2,3 (i.e., no change)
#'
#' @name ped_internal
NULL

#' @rdname ped_internal
#' @export
reorderPed = function(x, neworder = NULL, internal = FALSE) {
  if(is.pedList(x))
    stop2("Input is a ped list; reordering can only be done for a single component")
  if(!is.ped(x))
    stop2("Input is not a `ped` object or a list of such")
  if(is.singleton(x))
    return(x)

  if(is.null(neworder)) {
    orderby = if(hasNumLabs(x)) as.numeric(labels(x)) else labels(x)
    neworder = order(orderby)
    internal = TRUE
  }

  N = pedsize(x)
  if(anyDuplicated.default(neworder))
    stop2("Duplicated element of `neworder`: ", neworder(duplicated(neworder)))

  if(internal && !is.numeric(neworder))
    stop2("When `internal = TRUE`, `neworder` must be numeric: ", class(neworder)[1])

  # If not character, but not explicitly internal: interpret as labels if it makes sense
  if(!is.character(neworder) && !internal && all(neworder %in% x$ID))
    neworder = as.character(neworder)

  # Convert to internal indices
  if(is.character(neworder))
    neworder = internalID(x, neworder, errorIfUnknown = TRUE)

  ### By now: neworder assumed numeric with internal ordering

  if(!all(neworder %in% 1:N)) #
    stop2("Unknown index: ", .mysetdiff(neworder, 1:N))

  # If same order, return unchanged
  if(isTRUE(all.equal(neworder, 1:N)))
    return(x)

  # If subset of 1:N, expand:
  if(length(neworder) < N) {
    newo = neworder
    neworder = 1:N
    neworder[sort.int(newo, method = "quick")] = newo
  }

  # Convert to matrix with attributes
  xmatr = as.matrix(x)
  attr = attributes(xmatr)
  attr$LABELS = attr$LABELS[neworder]

  # Fix loop breakers
  if(!is.null(lp <- attr$LOOP_BREAKERS))
    attr$LOOP_BREAKERS = matrix(match(lp, neworder), ncol = 2)

  # Restore
  restorePed(xmatr[neworder, ], attrs = attr)
}


#' @rdname ped_internal
#' @export
parentsBeforeChildren = function(x) {
  if(is.pedList(x))
    return(lapply(x, parentsBeforeChildren))
  else if(!is.ped(x))
    stop2("Input is not a `ped` object or a list of such")
  if(is.singleton(x) || hasParentsBeforeChildren(x))
    return(x)

  neworder = 1:pedsize(x)
  i = 1
  while (i < pedsize(x)) {
    current = neworder[i]
    maxpar = max(match(c(x$FIDX[current], x$MIDX[current]), neworder, nomatch = 0))
    if (maxpar > i) { # push current indiv below below parents
      neworder[i:maxpar] = neworder[c((i+1):maxpar, i)]
    }
    else i = i + 1
  }
  reorderPed(x, neworder, internal = TRUE)
}

#' @rdname ped_internal
#' @export
hasParentsBeforeChildren = function(x) {
  if(is.pedList(x))
    return(all(vapply(x, hasParentsBeforeChildren, FUN.VALUE = logical(1))))
  else if(!is.ped(x))
    stop2("Input is not a `ped` object or a list of such")

  idx = 1:pedsize(x)
  father_before_child = x$FIDX < idx
  mother_before_child = x$MIDX < idx
  all(father_before_child & mother_before_child)
}


#' @rdname ped_internal
#' @export
foundersFirst = function(x) {
  if(is.pedList(x))
    return(lapply(x, foundersFirst))
  else if(!is.ped(x))
    stop2("Input is not a `ped` object or a list of such")

  fou = founders(x, internal = TRUE)

  # Check if all founders are already first
  if(length(fou) == max(fou))
    return(x)

  nonfou = nonfounders(x, internal = TRUE)
  reorderPed(x, neworder = c(fou, nonfou), internal = TRUE)
}




#' @rdname ped_internal
#' @export
internalID = function(x, ids, errorIfUnknown = TRUE) {
  if(is.pedList(x)) {
    comp = getComponent(x, ids, checkUnique = TRUE, errorIfUnknown = errorIfUnknown)
    idsInt = vapply(seq_along(ids), function(i) {
      if(is.na(comp[i])) NA_integer_ else internalID(x[[comp[i]]], ids[i])
    },
    FUN.VALUE = 1L)
    return(data.frame(id = ids, comp = comp, int = idsInt))
  }
  else if(!is.ped(x))
    stop2("Input is not a `ped` object or a list of such")

  idsInt = match(ids, labels(x))
  if (anyNA(idsInt) && errorIfUnknown)
    stop2("Unknown ID label: ", ids[is.na(idsInt)])
  idsInt
}

#' Get or modify pedigree labels
#'
#' Functions for getting or changing the ID labels of pedigree members.
#'
#' By default, `relabel(x)` relabels everyone as 1, 2, ..., in the order given
#' by the plot (top to bottom; left to right).
#'
#' Alternatively, `relabel(x, "generations")` labels the members in the top
#' generation I-1, I-2, ..., in the second generation II-1, II-2, ..., etc.
#'
#' @param x A `ped` object or a list of such.
#' @param new Either a character vector containing new labels, or one of the
#'   special words "asPlot" (default) or "generations". See Details and
#'   Examples.
#' @param old A vector of ID labels, of the same length as `new`. (Ignored if
#'   `new` is one of the special words.)
#' @param reorder A logical. If TRUE, [reorderPed()] is called on `x` after
#'   relabelling. Default: FALSE.
#' @param returnLabs A logical. If TRUE, the new labels are returned as a named
#'   character vector.
#' @param .alignment A list of alignment details for `x`, used if `new` equals
#'   "asPlot" or "generations". If not supplied, this is computed internally
#'   with [.pedAlignment()].
#'
#' @return
#'
#' * `labels()` returns a character vector containing the ID labels of all
#' pedigree members. If the input is a list of ped objects, the output is a list
#' of character vectors.
#'
#' * `relabel()` by default returns a `ped` object similar to `x`, but with
#'   modified labels. If `returnLabs` is TRUE, the new labels are returned as a
#'   named hcaracter vector
#'
#' @seealso [ped()]
#'
#' @examples
#'
#' x = nuclearPed()
#' x
#' labels(x)
#'
#' y = relabel(x, new = "girl", old = 3)
#' y
#'
#' # Back to the numeric labels
#' z = relabel(y)
#' stopifnot(identical(x,z))
#'
#' # Generation labels
#' relabel(x, "generations")
#'
#' @importFrom utils as.roman
#' @export
relabel = function(x, new = "asPlot", old = labels(x), reorder = FALSE, returnLabs = FALSE, .alignment = NULL) {
  if(is.list(old))
    old = unlist(old, use.names = FALSE)

  if(identical(new, "asPlot") || identical(new, "generations")) {

    if(is.pedList(x))
      stop2("`asPlot` cannot be used with ped lists")

    # Always reorder in this case
    reorder = TRUE

    p = .alignment$plist %||% .pedAlignment(x)$plist

    oldIdx = unlist(lapply(seq_along(p$n), function(i) p$nid[i, 1:p$n[i]]))

    # Remove duplicates
    dups = duplicated(oldIdx)
    if(any(dups))
      oldIdx = oldIdx[!dups]

    old = x$ID[oldIdx]

    if(identical(new, "generations")) {
      gen = rep(seq_along(p$n), p$n)

      if(any(dups))
        gen = gen[!dups]

      idx = unlist(lapply(1:max(gen), function(g) seq_len(sum(gen == g))))
      new = paste(as.roman(gen), idx, sep = "-")
    }
    else {
      new = as.character(seq_along(old))
    }

    names(new) = old
    if(returnLabs)
      return(new)
  }

  if(length(new) != length(old))
    stop2("Arguments `new` and `old` must have the same length")

  if(anyDuplicated.default(old) > 0)
    stop2("Duplicated entry in argument `old`: ", unique(old[duplicated(old)]))

  if(is.pedList(x)) {
    res = lapply(x, function(comp) {
      idx = old %in% labels(comp)
      if(!any(idx))
        comp
      else
        relabel(comp, new[idx], old[idx])
    })
    return(res)
  }
  else if(!is.ped(x))
    stop2("Input is not a `ped` object or a list of such")

  # Build new ID vector
  id = x$ID
  old_int = internalID(x, old)
  id[old_int] = new

  # Loop breakers
  if(!is.null(lb <- x$LOOP_BREAKERS)) {
    # Are any of the copies already changed by user?
    j = lb[, 'copy'] %in% old_int

    copy = lb[!j, 'copy']
    orig = lb[!j, 'orig']
    id[copy] = paste0("=", id[orig])
  }

  x$ID = id

  # Duplicated IDs after relabelling?
  if(anyDuplicated.default(id) > 0)
    stop2("Duplicated ID label: ", unique(id[duplicated(id)]))

  # Replace `pedmembers` attribute of each marker
  if(hasMarkers(x))
    x$MARKERS = lapply(x$MARKERS, `attr<-`, 'pedmembers', id)

  if(reorder)
    x = reorderPed(x)

  x
}

#' @param object A `ped` object
#' @param ... Not used
#'
#' @rdname relabel
#' @export
labels.ped = function(object, ...) {
  object$ID
}

#' @rdname relabel
#' @export
labels.list = function(object, ...) {
  if(is.pedList(object))
    lapply(object, labels.ped)
  else
    labels.default(object)
}

#' Get or set the sex of pedigree members
#'
#' Functions for retrieving or changing the sex of specified pedigree members.
#'
#' @param x A `ped` object or a list of such.
#' @param ids A character vector (or coercible to one) containing ID labels. If
#'   NULL, defaults to all members of `x`.
#' @param named A logical: return a named vector or not.
#' @param sex A numeric vector with entries 1 (= male), 2 (= female) or 0 (=
#'   unknown). If `ids` is NULL, `sex` must be named with ID labels. If `sex` is
#'   unnamed and shorter than `ids` it is recycled to `length(ids)`.
#' @param verbose A logical: Verbose output or not.
#'
#' @seealso [ped()]
#'
#' @return
#'
#' * `getSex(x, ids)` returns an integer vector of the same length as `ids`,
#' with entries 0 (unknown), 1 (male) or 2 (female).
#'
#' * `setSex(x, ids, sex)` returns a ped object similar to `x`, but where the
#' sex of `ids` is set according to the entries of `sex`
#'
#' * `swapSex(x, ids)` returns a ped object identical to `x`, but where the sex
#' of `ids` (and their spouses) are swapped (1 <-> 2).
#'
#' @examples
#' x = nuclearPed(father = "fa", mother = "mo", children = "ch")
#'
#' stopifnot(all.equal(
#'   getSex(x, named = TRUE),
#'   c(fa = 1, mo = 2, ch = 1)
#' ))
#'
#' # Make child female
#' setSex(x, ids = "ch", sex = 2)
#'
#' # Same, using a named vector
#' setSex(x, sex = c(ch = 2))
#'
#' # Swapping sex is sometimes easier,
#' # since spouses are dealt with automatically
#' swapSex(x, ids = "fa")
#'
#' # setting/getting sex in a pedlist
#' y = list(singleton(1, sex = 2), singleton(2), singleton(3))
#' sx = getSex(y, named = TRUE)
#' y2 = setSex(y, sex = sx)
#'
#' stopifnot(identical(y, y2))
#'
#' @importFrom stats setNames
#' @export
getSex = function(x, ids = NULL, named = FALSE) {
  if(is.pedList(x)) {
    sexVec = unlist(lapply(x, function(comp) comp$SEX), recursive = FALSE, use.names = FALSE)
    nms = unlist(labels(x), recursive = FALSE, use.names = FALSE)

    # Check for duplicates
    if(anyDuplicated.default(nms)) {
      dups = intersect(ids, nms[duplicated.default(nms)])
      if(length(dups))
        stop2("Duplicated ID label: ", dups)
    }

    names(sexVec) = nms
  }
  else
    sexVec = setNames(x$SEX, x$ID)

  storage.mode(sexVec) = "integer"

  res = if(is.null(ids)) sexVec else sexVec[as.character(ids)]

  if(!named)
    res = unname(res)

  res
}


#' @rdname getSex
#' @export
setSex = function(x, ids = NULL, sex) {
  if(!is.ped(x) && !is.pedList(x))
    stop2("Input is not a `ped` object or a list of such")

  if(is.null(ids))
    ids = names(sex)

  if(is.null(ids))
    stop2("If `ids` is NULL, then `sex` must be named")

  sex = as.integer(sex) # strips names, but ok since ids is already defined

  idsL = length(ids)
  sexL = length(sex)
  if(sexL < idsL)
    sex = rep(sex, length.out = idsL)
  else if(sexL > idsL)
    stop2("Argument `sex` is longer than `ids`")

  if(is.pedList(x)) {
    y = lapply(x, function(comp) {
      idsInt = internalID(comp, ids, errorIfUnknown = FALSE)
      notNA = !is.na(idsInt)
      comp$SEX[idsInt[notNA]] = sex[notNA]
      validatePed(comp)
      comp
    })
    return(y)
  }

  x$SEX[internalID(x, ids)] = sex

  validatePed(x)
  x
}



#' @rdname getSex
#' @export
swapSex = function(x, ids, verbose = TRUE) {

  if(is.pedList(x)) {
    y = lapply(x, function(comp) {
      idsComp = intersect(comp$ID, ids)
      if(!length(idsComp))
        comp
      else swapSex(comp, idsComp, verbose = verbose)
    })
    return(y)
  }
  else if(!is.ped(x))
    stop2("Input is not a `ped` object or a list of such")

  # Ignore individuals with unknown sex
  ids = ids[getSex(x, ids) != 0]
  if(!length(ids))
    return(x)

  idsInt = internalID(x, ids)
  FIDX = x$FIDX
  MIDX = x$MIDX
  spouses = c(MIDX[FIDX %in% idsInt], FIDX[MIDX %in% idsInt])

  if (!all(spouses %in% idsInt)) {
    extra = x$ID[setdiff(spouses, idsInt)]
    if (verbose)
      message("Changing sex of spouses as well: ", toString(extra))
    return(swapSex(x, union(ids, extra), verbose = verbose))
  }

  # Swap sex
  x$SEX[idsInt] = 3L - x$SEX[idsInt]

  # # Swap parents wherever any of the 'ids' occur as parents
  ids_as_parents = FIDX %in% idsInt # same with MIDX!
  FIDX[ids_as_parents] = x$MIDX[ids_as_parents]
  MIDX[ids_as_parents] = x$FIDX[ids_as_parents]
  x$FIDX = FIDX
  x$MIDX = MIDX

  x
}



#' Family identifier
#'
#' Functions for getting or setting the family ID of a `ped` object.
#'
#' @param x A `ped` object
#' @param value The new family ID, which must be (coercible to) a character
#'   string.
#' @param ... (Not used)
#'
#' @examples
#' x = nuclearPed(1)
#' famid(x) # empty string
#'
#' famid(x) = "trio"
#' famid(x)
#'
#' @export
`famid` = function(x, ...) {
  UseMethod("famid")
}

#' @rdname famid
#' @export
`famid.ped` = function(x, ...) {
  x$FAMID
}

#' @rdname famid
#' @export
`famid<-` = function(x, ..., value) {
  UseMethod("famid<-")
}

#' @rdname famid
#' @export
`famid<-.ped` = function(x, ..., value) {
  famid = as.character(value)
  if(length(famid) != 1) stop2("Replacement value must have length 1: ", famid)
  x$FAMID = famid
  x
}


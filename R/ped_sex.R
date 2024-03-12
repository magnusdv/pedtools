#' Get or set the sex of pedigree members
#'
#' Functions for retrieving or changing the sex of specified pedigree members.
#' When used in pedigree constructions, `swapSex()` is usually more convenient
#' than `setSex()`, since it deals with spouses automatically.
#'
#' To set unknown sex, use `setSex(x, ids, sex = 0)`. Note that if a nonfounder
#' has unknown sex the pedigree cannot be plotted in the usual way, only with
#' `plot(x, arrows = TRUE)`.
#'
#' @param x A `ped` object or a list of such.
#' @param ids A vector identifying members of `x`, or a function, in which case
#'   it is replaced with `ids(x)` labels. If NULL, defaults to all members of
#'   `x`.
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
#' of `ids` (and their spouses) are swapped (1 <-> 2). Individuals of unknown
#' sex are ignored.
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
#' # Same, using a function (setting all leaves to be female)
#' setSex(x, ids = leaves, sex = 2)
#'
#' # swapSex() deals with spouses automatically
#' swapSex(x, ids = "fa")
#'
#' # setting/getting sex in a pedlist
#' y = singletons(id = 1:3, sex = c(2,1,1))
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
    nms = labels(x)

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

  ispedlist <- is.pedList(x)
  if(!is.ped(x) && !ispedlist)
    stop2("Input is not a `ped` object or a list of such")

  if(is.function(ids))
    ids = ids(x)
  else
    ids = ids %||% names(sex) %||% stop2("If `ids` is NULL, then `sex` must be named")

  sex = as.integer(sex) # strip names

  idsL = length(ids)
  sexL = length(sex)
  if(sexL < idsL)
    sex = rep(sex, length.out = idsL)
  else if(sexL > idsL)
    stop2("Argument `sex` is longer than `ids`")

  if(ispedlist) {
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

  ispedlist <- is.pedList(x)
  if(!is.ped(x) && !ispedlist)
    stop2("Input is not a `ped` object or a list of such")

  if(is.function(ids))
    ids = ids(x)

  if(ispedlist) {
    y = lapply(x, function(comp) {
      idsComp = intersect(comp$ID, ids)
      if(!length(idsComp))
        comp
      else swapSex(comp, idsComp, verbose = verbose)
    })
    return(y)
  }

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

  # Update 'sex' attribute of each marker
  for(i in seq_along(x$MARKERS))
    attr(x$MARKERS[[i]], "sex") = x$SEX

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


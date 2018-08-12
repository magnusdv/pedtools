#' Pedigree modifications
#'
#' Functions for modifying various attributes of a 'ped' object.
#'
#' @param x A `ped` object.
#' @param ids A character (or coercible to character) with ID labels of one or
#'   more pedigree members.
#' @param labels A character (or coercible to character) of length `pedsize(x)`
#' @param new,old Character vectors (or coercible to character) of the same
#'   length. ID labels in `old` are replaced by those in `new`.
#' @param verbose A logical: Verbose output or not.
#'
#' @return The modified `ped` object.
#' @author Magnus Dehli Vigeland
#' @seealso [ped()], [ped_add]
#'
#' @examples
#'
#' x = nuclearPed(1)
#'
#' # To see the effect of each command below, use plot(x) in between.
#' x = swapSex(x, 3)
#' x = relabel(x, new="girl", old=3)
#'
#' @name ped_modify
NULL

#' @rdname ped_modify
#' @export
swapSex = function(x, ids, verbose = TRUE) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(!length(ids)) return(x)
  ids = internalID(x, ids)
  labs = labels(x)
  FID = x$FID
  MID = x$MID
  spouses = c(MID[FID %in% ids], FID[MID %in% ids])

  if (!all(spouses %in% ids)) {
    if (verbose) {
      extra = setdiff(spouses, ids)
      message("Changing sex of spouses as well: ", toString(labs[extra]))
    }
    return(swapSex(x, labs[union(ids, spouses)]))
  }

  # Swap sex
  x$SEX[ids] = 3 - x$SEX[ids]

  # # Swap parents wherever any of the 'ids' occur as parents
  ids_as_parents = FID %in% ids # same with MID!
  FID[ids_as_parents] = x$MID[ids_as_parents]
  MID[ids_as_parents] = x$FID[ids_as_parents]
  x$FID = FID
  x$MID = MID

  x
}

#' @rdname ped_modify
#' @export
relabel = function(x, new, old=labels(x)) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(length(new) != length(old))
    stop2("Arguments `new` and `old` must have the same length")
  xlabs = labels(x)
  old_idx = match(old, xlabs)
  if(anyNA(old_idx))
    stop2("Unknown ID label: ", old[is.na(old_idx)])

  x$LABELS[old_idx] = new

  if(hasMarkers(x)) {
    x$markerdata = lapply(x$markerdata, `attr<-`, 'pedmembers', x$LABELS)
  }
  x
}


#' @rdname ped_modify
#' @export
setLabels = function(x, labels) {
  message("`setLabels()` is deprecated. Use `relabel()` instead")
  labels = as.character(labels)
  #assert_that(is.ped(x), length(labels) == pedsize(x), !anyDuplicated(labels))
  x$LABELS = labels
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


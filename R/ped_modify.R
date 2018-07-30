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
#' @param famid A character of length 1. If missing, an emtpy string is used.
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
#' x = setFamid(x, "Family 1")
#'
#' @name ped_modify
NULL

#' @rdname ped_modify
#' @export
swapSex = function(x, ids, verbose = TRUE) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(!length(ids)) return(x)
  ids = internalID(x, ids)
  FID = x$FID
  MID = x$MID
  spouses = c(MID[FID %in% ids], FID[MID %in% ids])

  if (!all(spouses %in% ids)) {
    if (verbose) {
      extra = setdiff(spouses, ids)
      message("Changing sex of spouses as well: ", toString(x$LABELS[extra]))
    }
    return(swapSex(x, x$LABELS[union(ids, spouses)]))
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
relabel = function(x, new, old=x$LABELS) {
  assert_that(is.ped(x), all(old %in% x$LABELS), length(new)==length(old))
  lab = x$LABELS
  lab[match(old, lab)] = new
  x$LABELS = lab

  if(hasMarkers(x)) {
    x$markerdata[] = lapply(x$markerdata, `attr<-`, 'pedmembers', lab)
  }
  x
}


#' @rdname ped_modify
#' @export
setLabels = function(x, labels) {
  labels = as.character(labels)
  assert_that(is.ped(x), length(labels) == pedsize(x), !anyDuplicated(labels))
  x$LABELS = labels
  x
}


#' @rdname ped_modify
#' @export
setFamid = function(x, famid) {
  famid = as.character(famid)
  if(!length(famid)) famid=""
  assert_that(is.ped(x), length(famid) == 1)
  x$FAMID = famid
  x
}


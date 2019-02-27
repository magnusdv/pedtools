#' Inbreeding coefficients of founders
#'
#' Functions to get or set inbreeding coefficients for the pedigree founders.
#'
#' @param x A `ped` object
#' @param ids Any subset of `founders(x)`. If `ids` is missing in
#'   `founderInbreeding()`, it is set to `founders(x)`.
#' @param named A logical: If TRUE, the output vector is named with the ID
#'   labels
#'
#' @return For `founderInbreeding`, a numeric vector of the same length as
#' `ids`, containing the founder inbreeding coefficients.
#'
#' For `founderInbreeding<-` the updated `ped` object is returned.
#'
#' @examples
#' x = nuclearPed(1)
#' founderInbreeding(x, ids = '1') = 1
#' founderInbreeding(x, named = TRUE)
#'
#' # Setting all founders at once (replacement value is recycled)
#' founderInbreeding(x, ids = founders(x)) = 0.5
#' founderInbreeding(x, named = TRUE)
#'
#' # Alternative syntax, using a named vector
#' founderInbreeding(x) = c('1'=0.1, '2'=0.2)
#' founderInbreeding(x, named = TRUE)
#'
#' @export
founderInbreeding = function(x, ids, named = FALSE) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  fou = founders(x)

  if (missing(ids))
    ids = fou
  else if(any(!ids %in% fou)) {
    internalID(x, ids) # quick hack to catch unknown labels
    stop2("Pedigree member is not a founder: ", setdiff(ids, fou))
  }

  finb = x$FOUNDER_INBREEDING

  if(is.null(finb))
    finb = rep(0, length(fou))

  res = finb[match(ids, fou)]
  if(named)
    names(res) = ids
  res
}

#' @param value A numeric of the same length as `ids`, entries in the interval
#'   `[0, 1]`. If the vector is named, then the names are interpreted as ID
#'   labels of the founders whose inbreeding coefficients should be set. In this
#'   case, the `ids` argument should not be used. (See examples.)
#'
#'
#' @rdname founderInbreeding
#' @export
`founderInbreeding<-` = function(x, ids, value) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(!is.numeric(value))
    stop2("Inbreeding coefficients must be numeric: ", value)

  illegal = value < 0 | value > 1
  if(any(illegal))
    stop2("Inbreeding coefficients must be in the interval [0, 1]: ", value[illegal])

  if(missing(ids) && is.null(names(value)))
    stop2("When argument `ids` is missing, the replacement vector must be named")
  if(!missing(ids) && !is.null(names(value)))
    stop2("When the replacement vector is named, the `ids` argument must be missing")

  if(missing(ids)) ids = names(value)

  if(length(value) == 1)
    value = rep_len(value, length(ids))
  else if(length(ids) != length(value))
    stop2("Replacement vector must have length 1 or length(ids)")

  if(anyDuplicated.default(ids) > 0)
    stop2("Duplicated ID label: ", ids[duplicated(ids)])

  fou = founders(x)
  if(any(!ids %in% fou)) {
    internalID(x, ids) # quick hack to catch unknown labels
    stop2("Pedigree member is not a founder: ", setdiff(ids, fou))
  }

  finb = x$FOUNDER_INBREEDING
  if(is.null(finb))
    finb = rep(0, length(fou))

  finb[match(ids, fou)] = value
  x$FOUNDER_INBREEDING = finb
  x
}

#' Inbreeding coefficients of founders
#'
#' Functions to get or set inbreeding coefficients for the pedigree founders.
#'
#' @param x A `ped` object.
#' @param ids Any subset of `founders(x)`. If `ids` is missing in
#'   `founderInbreeding()`, it is set to `founders(x)`.
#' @param named A logical: If TRUE, the output vector is named with the ID
#'   labels.
#' @param chromType Either "autosomal" (default) or "x".
#'
#' @return For `founderInbreeding`, a numeric vector of the same length as
#'   `ids`, containing the founder inbreeding coefficients.
#'
#'   For `founderInbreeding<-` the updated `ped` object is returned.
#'
#' @examples
#' x = nuclearPed(father = "fa", mother = "mo", child = 1)
#' founderInbreeding(x, "fa") = 1
#' founderInbreeding(x, named = TRUE)
#'
#' # Setting all founders at once (replacement value is recycled)
#' founderInbreeding(x, ids = founders(x)) = 0.5
#' founderInbreeding(x, named = TRUE)
#'
#' # Alternative syntax, using a named vector
#' founderInbreeding(x) = c(fa = 0.1, mo = 0.2)
#' founderInbreeding(x, named = TRUE)
#'
#' @export
founderInbreeding = function(x, ids, named = FALSE, chromType = "autosomal") {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(!chromType %in% c("autosomal", "x"))
    stop2("Argument `chromType` must be a either 'autosomal' or 'x': ", chromType)

  finb = x$FOUNDER_INBREEDING[[chromType]]

  if(is.null(finb)) {
    isFou = x$FIDX == 0
    finb = rep_len(0, sum(isFou))
    if(chromType == "x") {
      maleFou = x$SEX[isFou] == 1
      finb[maleFou] = 1 # always 1 for males
    }
  }

  # Quick return if `ids` is missing (so no ID checks are needed)
  if(missing(ids)) {
    if(named)
      names(finb) = founders(x)
    return(finb)
  }

  fou = founders(x)
  if(any(!ids %in% fou)) {
    internalID(x, ids) # quick hack to catch unknown labels
    stop2("Pedigree member is not a founder: ", setdiff(ids, fou))
  }

  finb = finb[match(ids, fou)]
  if(named)
    names(finb) = ids

  finb
}

#' @param value A numeric of the same length as `ids`, entries in the interval
#'   `[0, 1]`. If the vector is named, then the names are interpreted as ID
#'   labels of the founders whose inbreeding coefficients should be set. In this
#'   case, the `ids` argument should not be used. (See examples.)
#'
#'
#' @rdname founderInbreeding
#' @export
`founderInbreeding<-` = function(x, ids, chromType = "autosomal", value) {
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

  chromType = match.arg(tolower(chromType), c("autosomal", "x"))
  if(chromType =="x" && any(value[ids %in% males(x)] < 1))
    stop2("X chromosomal inbreeding coefficient of males cannot be less than 1: ", value[ids %in% males(x)])

  # Back compatibility: Allow previous version where FOUNDER_INBREEDING was just a vector
  finb = x$FOUNDER_INBREEDING
  if(!is.list(finb) && is.numeric(finb))
    finb = list(autosomal = finb, x = NULL)

  current = finb[[chromType]]

  if(is.null(current))
    current = switch(chromType,
                     autosomal = rep(0, length(fou)),
                     x = ifelse(getSex(x, fou) == 1, 1, 0))

  current[match(ids, fou)] = value
  x$FOUNDER_INBREEDING[[chromType]] = current
  x
}

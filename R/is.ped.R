#' Is an object a `ped` object?
#'
#' Functions for checking whether an object is a [ped()] object, a [singleton()]
#' or a list of such.
#'
#' Note that the `singleton` class inherits from `ped`, so if `x` is a
#' singleton, `is.ped(x)` returns TRUE.
#'
#' @param x Any `R` object.
#' @return For `is.ped()`: TRUE if `x` is a `ped` or `singleton` object, otherwise FALSE.
#'
#'   For `is.singleton()`: TRUE if `x` is a `singleton` object, otherwise FALSE.
#'
#'   For `is.pedList()`: TRUE if `x` is a list of `ped` and/or `singleton`
#'   objects, otherwise FALSE.
#'
#' @author Magnus Dehli Vigeland
#' @seealso [ped()]
#' @examples
#'
#' x1 = nuclearPed(1)
#' x2 = singleton(1)
#' stopifnot(is.ped(x1), !is.singleton(x1),
#'           is.ped(x2), is.singleton(x2),
#'           is.pedList(list(x1,x2)))
#'
#' @export
is.ped = function(x)
  inherits(x, "ped")

#' @rdname is.ped
#' @export
is.singleton = function(x)
  inherits(x, "singleton")

#' @rdname is.ped
#' @export
is.pedList = function(x) {
  # Most common FALSE case: ped
  if(is.ped(x) || !is.list(x) || length(x) == 0)
    return(FALSE)

  if(inherits(x, "pedList"))
    return(TRUE)

  for(comp in x) {
    if(!inherits(comp, what = "ped"))
      return(FALSE)
  }
  return(TRUE)
}

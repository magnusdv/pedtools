#' Is an object a ped object?
#'
#' Functions for checking whether an object is a [ped()] object, a [singleton()] or
#' a list of such.
#'
#' Note that the `singleton` class inherits from `ped`, so if
#' `x` is a singleton, `is.ped(x)` returns TRUE.
#'
#' @param x Any `R` object.
#' @return For `is.ped`: TRUE if `x` is a ped (or singleton)
#' object, and FALSE otherwise.\cr For `is.singleton`: TRUE if `x` is
#' a singleton object, and FALSE otherwise.\cr For `is.ped.list`: TRUE
#' if `x` is a list of ped/singleton objects.
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
    isTRUE(is.list(x) &&
           length(x) > 0 &&
           all(vapply(x, function(comp) inherits(comp, what = "ped"), FUN.VALUE = logical(1))))
}

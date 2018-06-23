#' Is an object a ped object?
#' 
#' Functions for checking whether an object is a \code{\link{ped}} object, a \code{\link{singleton}} or
#' a list of such.
#' 
#' Note that the \code{singleton} class inherits from \code{ped}, so if
#' \code{x} is a singleton, \code{is.ped(x)} returns TRUE.
#' 
#' @param x Any \code{R} object.
#' @return For \code{is.ped}: TRUE if \code{x} is a ped (or singleton)
#' object, and FALSE otherwise.\cr For \code{is.singleton}: TRUE if \code{x} is
#' a singleton object, and FALSE otherwise.\cr For \code{is.ped.list}: TRUE
#' if \code{x} is a list of ped/singleton objects.
#' @author Magnus Dehli Vigeland
#' @seealso \code{\link{ped}}
#' @examples
#' 
#' x1 = nuclearPed(1)
#' x2 = singleton(1)
#' stopifnot(is.ped(x1), !is.singleton(x1), 
#'           is.ped(x2), is.singleton(x2),
#'           is.ped.list(list(x1,x2)))
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
is.ped.list = function(x) 
    isTRUE(is.list(x) && all(sapply(x, inherits, "ped")))

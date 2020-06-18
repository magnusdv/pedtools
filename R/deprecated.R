#' #' Deprecated functions
#' #'
#' #' These functions have been renamed (from "snake_case" to "camelCase") and will eventually be removed.
#' #'
#' #' @param x object
#' #'
#' #' @name deprecated
#' NULL
#'
#' #' @rdname deprecated
#' #' @export
#' has_common_ancestor = function(x) {
#'   warning("The function `has_common_ancestor()` is renamed to `hasCommonAncestor()` and will eventually be removed",
#'           call. = FALSE)
#'   hasCommonAncestor(x)
#' }

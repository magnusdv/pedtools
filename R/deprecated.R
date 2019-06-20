#' Deprecated functions
#'
#' These functions have been renamed (from "snake_case" to "camelCase") and will eventually be removed.
#'
#' @param x object
#'
#' @name deprecated
NULL

#' @rdname deprecated
#' @export
has_common_ancestor = function(x) {
  warning("The function `has_common_ancestor()` is renamed to `hasCommonAncestor()` and will eventually be removed",
          call. = FALSE)
  hasCommonAncestor(x)
}

#' @rdname deprecated
#' @export
has_unbroken_loops = function(x) {
  warning("The function `has_unbroken_loops()` is renamed to `hasUnbrokenLoops()` and will eventually be removed",
          call. = FALSE)
  hasUnbrokenLoops(x)
}

#' @rdname deprecated
#' @export
has_selfing = function(x) {
  warning("The function `has_selfing()` is renamed to `hasSelfing()` and will eventually be removed",
          call. = FALSE)
  hasSelfing(x)
}

#' @rdname deprecated
#' @export
has_inbred_founders = function(x) {
  warning("The function `has_inbred_founders()` is renamed to `hasInbredFounders()` and will eventually be removed",
          call. = FALSE)
  hasInbredFounders(x)
}

#' @rdname deprecated
#' @export
has_parents_before_children = function(x) {
  warning("The function `has_parents_before_children()` is renamed to `hasParentsBeforeChildren()` and will eventually be removed",
          call. = FALSE)
  hasParentsBeforeChildren(x)
}

#' @rdname deprecated
#' @export
parents_before_children = function(x) {
  warning("The function `parents_before_children()` is renamed to `parentsBeforeChildren()` and will eventually be removed",
          call. = FALSE)
  parentsBeforeChildren(x)
}

#' @rdname deprecated
#' @export
is_Xmarker = function(x) {
  warning("The function `is_Xmarker()` is renamed to `isXmarker()` and will eventually be removed",
          call. = FALSE)
  isXmarker(x)
}

#' @rdname deprecated
#' @export
restore_ped = function(x) {
  warning("The function `restore_ped()` is renamed to `restorePed()` and will eventually be removed",
          call. = FALSE)
  restorePed(x)
}

#' @rdname deprecated
#' @export
validate_ped = function(x) {
  warning("The function `validate_ped()` is renamed to `validatePed()` and will eventually be removed",
          call. = FALSE)
  validatePed(x)
}

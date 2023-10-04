#' Get or modify pedigree labels
#'
#' Functions for getting or changing the ID labels of pedigree members.
#'
#' By default, `relabel(x)` relabels everyone as 1, 2, ..., in the order given
#' by the plot (top to bottom; left to right).
#'
#' Alternatively, `relabel(x, "generations")` labels the members in the top
#' generation I-1, I-2, ..., in the second generation II-1, II-2, ..., etc.
#'
#' @param x A `ped` object or a list of such.
#' @param new The following values are valid (see Details and Examples):
#'   * a character vector containing new labels. If named, interpreted as
#'   `old = new`
#'   * a function, which should take the old labels as input and output a
#'   character of the same length
#'   * one of the special keywords "asPlot" (default) or "generations"
#' @param old A vector of ID labels, of the same length as `new`. (Ignored if
#'   `new` is one of the special words.) If not given, taken from the names of
#'   `new` if these exist.
#' @param reorder A logical. If TRUE, [reorderPed()] is called on `x` after
#'   relabelling. Default: FALSE.
#' @param returnLabs A logical. If TRUE, the new labels are returned as a named
#'   character vector.
#' @param .alignment A list of alignment details for `x`, used if `new` equals
#'   "asPlot" or "generations". If not supplied, this is computed internally
#'   with [.pedAlignment()].
#'
#' @return
#'
#' * `labels()` returns a character vector containing the ID labels of all
#' pedigree members. If the input is a list of ped objects, the output is a list
#' of character vectors.
#'
#' * `relabel()` by default returns a `ped` object similar to `x`, but with
#' modified labels. If `returnLabs` is TRUE, the new labels are returned as a
#' named character vector
#'
#' @seealso [ped()]
#'
#' @examples
#'
#' x = nuclearPed()
#' x
#' labels(x)
#'
#' y = relabel(x, new = "girl", old = 3)
#' y
#'
#' # Back to the numeric labels
#' z = relabel(y)
#' stopifnot(identical(x,z))
#'
#' # Generation labels
#' relabel(x, "generations")
#'
#' @importFrom utils as.roman
#' @export
relabel = function(x, new = "asPlot", old = labels(x), reorder = FALSE,
                   returnLabs = FALSE, .alignment = NULL) {
  if(is.list(old))
    old = unlist(old, use.names = FALSE)

  if(identical(new, "asPlot") || identical(new, "generations")) {

    if(is.pedList(x))
      stop2("`asPlot` cannot be used with ped lists")

    # Always reorder in this case
    reorder = TRUE

    p = .alignment$plist %||% .pedAlignment(x)$plist

    oldIdx = unlist(lapply(seq_along(p$n), function(i) p$nid[i, 1:p$n[i]]))

    # Remove duplicates
    dups = duplicated(oldIdx)
    if(any(dups))
      oldIdx = oldIdx[!dups]

    # Check for misalignment
    if(length(oldIdx) != length(x$ID)) {
      warning("Alignment error; cannot relabel this pedigree according to plot order")
      return(x)
    }

    old = x$ID[oldIdx]

    if(identical(new, "generations")) {
      gen = rep(seq_along(p$n), p$n)

      if(any(dups))
        gen = gen[!dups]

      idx = unlist(lapply(1:max(gen), function(g) seq_len(sum(gen == g))))
      new = paste(as.roman(gen), idx, sep = "-")
    }
    else {
      new = as.character(seq_along(old))
    }

    names(new) = old
    if(returnLabs)
      return(new)
  }

  if(missing(old) && !is.null(names(new)))
    old = names(new)

  if(is.function(new))
    new = new(old)

  if(length(new) != length(old))
    stop2("Arguments `new` and `old` must have the same length")

  if(anyDuplicated.default(old) > 0)
    stop2("Duplicated entry in argument `old`: ", unique(old[duplicated(old)]))

  if(is.pedList(x)) {
    res = lapply(x, function(comp) {
      idx = old %in% labels(comp)
      if(!any(idx))
        comp
      else
        relabel(comp, new[idx], old[idx])
    })
    return(res)
  }
  else if(!is.ped(x))
    stop2("Input is not a `ped` object or a list of such")

  # Build new ID vector
  id = x$ID
  old_int = internalID(x, old)
  id[old_int] = new

  # Loop breakers
  if(!is.null(lb <- x$LOOP_BREAKERS)) {
    # Are any of the copies already changed by user?
    j = lb[, 'copy'] %in% old_int

    copy = lb[!j, 'copy']
    orig = lb[!j, 'orig']
    id[copy] = paste0("=", id[orig])
  }

  x$ID = id

  # Duplicated IDs after relabelling?
  if(anyDuplicated.default(id) > 0)
    stop2("Duplicated ID label: ", unique(id[duplicated(id)]))

  # Replace `pedmembers` attribute of each marker
  if(hasMarkers(x))
    x$MARKERS = lapply(x$MARKERS, `attr<-`, 'pedmembers', id)

  if(reorder)
    x = reorderPed(x)

  x
}

#' @param object A `ped` object
#' @param ... Not used
#'
#' @rdname relabel
#' @export
labels.ped = function(object, ...) {
  object$ID
}

#' @rdname relabel
#' @export
labels.list = function(object, ...) {
  if(is.pedList(object))
    lapply(object, labels.ped)
  else
    labels.default(object)
}

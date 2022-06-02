#' Merge two pedigrees
#'
#' This function merges two `ped` objects, joining them at the indicated
#' individuals. Only ped objects without marker data are supported.
#'
#' Some internal checks are done to ensure that merging individuals have the
#' same sex and the same parents.
#'
#' If `relabel = FALSE`, some relabelling might still be performed in order to
#' ensure unique labels for everyone. Specifically, this is the case if some ID
#' labels occur in both `x` and `y` other than those given in the `by` argument.
#' In such cases, the relevant members of `y` get a suffix `.y`.
#'
#' @param x,y [ped()] objects
#' @param by The individuals to merge by. The most general form uses a named
#'   vector with entries of the form `id.x = id.y` (see Examples). If the vector
#'   is unnamed, it is assumed that the merging individuals have the same labels
#'   in both pedigrees. Finally, if `by = NULL` (default), it is set to
#'   `intersect(labels(x), labels(y))`.
#' @param relabel A logical, by default FALSE. If TRUE, `relabel(..., "asPlot")`
#'   is run on the merged pedigree before returning.
#' @param ... further arguments passed along to [ped()], e.g. `famid`,
#'   `validate` and `reorder`.
#'
#' @return A `ped` object.
#' @author Magnus Dehli Vigeland
#'
#' @examples
#'
#' ############
#' # Example 1
#' # A family trio where each parent have first cousin parents.
#' ############
#'
#' # Trio
#' x = nuclearPed(1)
#'
#' # Add paternal family
#' pat = cousinPed(1, child = TRUE)
#' x = mergePed(x, pat, by = c("1" = "9"))
#'
#' # Maternal family
#' mat = cousinPed(1, child = TRUE) |> swapSex("9")
#' x = mergePed(x, mat, by = c("2" = "9"))
#'
#' # Relabel (Alternative: add `relabel = TRUE` in the previous call)
#' x = relabel(x, "asPlot")
#'
#' plot(x)
#'
#'
#' ##################################
#' # Example 2: Double first cousins
#' ##################################
#'
#' # First cousins, whose fathers are brothers
#' y = cousinPed(degree = 1)
#'
#' # Create two sisters
#' motherPed = nuclearPed(2, sex = 2)
#'
#' # Plot to see who is who: `plotPedList(list(y, motherPed))`
#'
#' # Merge
#' z = mergePed(y, motherPed, by = c("4" = 3, "6" = 4), relabel = TRUE)
#' plot(z)
#'
#' @export
mergePed = function(x, y, by = NULL, relabel = FALSE, ...) {
  if (!is.null(x$MARKERS) || !is.null(y$MARKERS))
    stop2("Merging is only supported for pedigrees without marker data")

  xlabs = labels(x)
  ylabs = labels(y)

  if(is.null(by))
    by.x = by.y = intersect(xlabs, ylabs)
  else {
    by.x = names(by) %||% as.character(by)
    by.y = as.character(by)
  }

  if(!length(by.x))
    stop2("Please indicate merges with the `by` argument: `by = c(<idx> = <idy>, ...)`")

  # Check IDs
  if(!all(as.character(by.x) %in% xlabs))
    stop2("Unknown ID label in pedigree 1: ", setdiff(by.x, xlabs))
  if(!all(as.character(by.y) %in% ylabs))
    stop2("Unknown ID label in pedigree 2: ", setdiff(by.y, ylabs))

  # Check genders
  sameSex = getSex(x, by.x) == getSex(y, by.y)
  if(!all(sameSex)) {
    mess = sapply(which(!sameSex), function(i) {
      id1 = by.x[i]; id2 = by.y[i]
      sex1 = c("male", "female")[getSex(x, id1)]
      sex2 = c("male", "female")[getSex(y, id2)]
      sprintf("\n  '%s' in pedigree 1 (%s)  <=/=>  '%s' in pedigree 2 (%s)", id1, sex1, id2, sex2)
    })
    stop2("Gender mismatch", mess)
  }

  # Non-merging members of y that need new labels
  dups = intersect(xlabs, setdiff(y$ID, by.y))
  while(length(dups)) {
    y = relabel(y, old = dups, new = paste0(dups, ".y"))
    dups = intersect(xlabs, setdiff(y$ID, by.y))
  }

  # Relabel 'by.y'
  y = relabel(y, old = by.y, new = by.x)

  # Re-extract labs
  ylabs = labels(y)

  del = list(x = numeric(), y = numeric())
  for (i in by.x) {
    parx = parents(x, i)
    pary = parents(y, i)

    if (length(pary) == 0)
      del$y = c(del$y, i)
    else if (length(parx) == 0)
      del$x = c(del$x, i)
    else if (all(parx == pary))
      del$y = c(del$y, i)
    else {
      mess = sprintf("Parent mismatch for individual %s.\nDid you forget to include some pairs in the `by` argument?", i)
      stop2(mess)
    }
  }

  xm = as.data.frame(x)
  ym = as.data.frame(y)
  zm = rbind(xm[!xlabs %in% del$x, ], ym[!ylabs %in% del$y, ])

  z = ped(id = zm$id, fid = zm$fid, mid = zm$mid, sex = zm$sex, ...)
  if(relabel)
    z = relabel(z, "asPlot")

  z
}

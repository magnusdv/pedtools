#' Connected pedigree components
#'
#' Compute the connected parts of a pedigree. This is an important step when
#' converting pedigree data from other formats (where disconnected pedigrees may
#' be allowed) to `pedtools` (which requires pedigrees to be connected).
#'
#' @param id A vector of ID labels (character or numeric).
#' @param fid The ID labels of the fathers (or "0" if missing).
#' @param mid The ID labels of the mothers (or "0" if missing).
#' @param fidx,midx (For internal use mostly.) Integer vectors with paternal
#'   (resp. maternal) indices. These may be given instead of `id`, `fid`, `mid`.
#' @return A list, where each element is a subset of `id` constituting a
#'   connected pedigree.
#'
#' @examples
#'
#' # A trio (1-3) and a singleton (4)
#' x = data.frame(id = 1:4, fid = c(2,0,0,0), mid = c(3,0,0,0))
#' connectedComponents(x$id, x$fid, x$mid)
#'
#' @export
connectedComponents = function(id, fid = NULL, mid = NULL, fidx = NULL, midx = NULL) {
  seqn = seq_along(id)

  if(is.null(fidx)) {
    fidx = match(fid, id, nomatch = 0L)
    midx = match(mid, id, nomatch = 0L)
  }

  n = length(id)

  # Fast path: all singletons
  if(!any(fidx) && !any(midx))
    return(unname(as.list(id)))

  getComp = function(start) {
    seen = logical(n)
    seen[start] = TRUE

    for(i in seqn) {
      old = seen

      # Upwards: add parents of seen individuals
      seen[fidx[old]] = TRUE
      seen[midx[old]] = TRUE

      # Downwards: add children with at least one seen parent
      parSeen = c(FALSE, seen)
      seen = seen | parSeen[fidx + 1L] | parSeen[midx + 1L]

      if(all(seen) || identical(seen, old))
        break
    }

    seen
  }

  # Fast path: single connected component
  seen = getComp(1L)
  if(all(seen))
    return(list(id))

  # Fast path: one component plus remaining singletons
  rest = !seen
  if(!any(fidx[rest]) && !any(midx[rest]))
    return(c(list(id[seen]), unname(as.list(id[rest]))))

  comp = integer(n)
  tag = 1L
  comp[seen] = tag

  repeat {
    start = match(0L, comp)
    if(is.na(start))
      break

    tag = tag + 1L
    comp[getComp(start)] = tag
  }

  unname(split.default(id, comp))
}

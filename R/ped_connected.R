#' Connected pedigree components
#'
#' Compute the connected parts of a pedigree. This is an important step when
#' converting pedigree data from other formats (where disconnected pedigrees may
#' be allowed) to `pedtools` (which requires pedigrees to be connected).
#'
#' @param id A vector of ID labels (character or numeric)
#' @param fid The ID labels of the fathers (or "0" if missing)
#' @param mid The ID labels of the mothers (or "0" if missing)
#'
#' @return A list, where each element is a subset of `id` constituting a
#'   connected pedigree
#'
#' @examples
#' # A trio (1-3) and a singleton (4)
#' x = data.frame(id = 1:4, fid = c(2,0,0,0), mid = c(3,0,0,0))
#' connectedComponents(x$id, x$fid, x$mid)
#'
#' @export
connectedComponents = function(id, fid, mid) {
  # Placeholder for final components
  comps = list()

  # Starting point: List of all founders and trios
  temp = lapply(seq_along(id), function(i) .mysetdiff(c(id[i], fid[i], mid[i]), 0))

  while (length(temp) > 0) {
    # Check if first vector overlaps with any of the others
    a = temp[[1]]
    remov = numeric()
    for (j in seq_along(temp)[-1]) {
      if (any(match(a, temp[[j]], nomatch = 0) > 0)) {
        a = unique.default(c(a, temp[[j]]))
        remov = c(remov, j)
      }
    }

    if (length(remov) > 0) {
      # Remove any overlapping vectors, and insert the union as first element
      temp[remov] = NULL
      temp[[1]] = a
    } else {
      # If no overlaps, we have a maximal component. Move to comps and remove from temp.
      comps = c(comps, list(sort.default(a)))
      temp[[1]] = NULL
    }
  }
  comps
}

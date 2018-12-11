#' Connected pedigree components
#'
#' Compute the connected parts of a pedigree. This is an important step when
#' converting pedigree data from other formats (where disconnected pedigrees may
#' be allowed) to `pedtools` (which requires pedigrees to be connected).
#'
#' @param id A vector of ID labels (character or numeric)
#' @param fid The ID labels of the fathers (or "0" if missing)
#' @param mid The ID labels of the mothers (or "0" if missing)
#' @param fidx,midx (For internal use mostly) Integer vectors with paternal
#'   (resp maternal) indices. These may be given instead of `id`, `fid`, `mid`.
#' @return A list, where each element is a subset of `id` constituting a
#'   connected pedigree
#'
#' @examples
#' # A trio (1-3) and a singleton (4)
#' x = data.frame(id = 1:4, fid = c(2,0,0,0), mid = c(3,0,0,0))
#' connectedComponents(x$id, x$fid, x$mid)
#'
#' @export
connectedComponents = function(id, fid, mid, fidx=NULL, midx=NULL) {
  if(!missing(id)) {
    fidx = match(fid, id, nomatch = 0L)
    midx = match(mid, id, nomatch = 0L)
  }
  else {
    id = seq_along(fidx)
  }
  seqn = seq_along(id)

  adjacencyList = lapply(seqn, function(i) {
    fa = fidx[i]; mo = midx[i]
    c(if(fa > 0) fa, if(mo > 0 && mo != fa) mo, seqn[fidx == i | midx == i])
  })

  env = new.env()
  env$comp = integer(length(fidx))

  DFS = function(i, tag) {
    env$comp[i] = tag
    for(j in adjacencyList[[i]]) {
      if(env$comp[j] == 0)
        DFS(j, tag)
    }
  }

  founders = seqn[fidx == 0 & midx == 0]
  tag = 0
  for(fou in founders) {
    if(env$comp[fou] == 0) {
      tag = tag + 1
      DFS(fou, tag)
    }
  }

  comps = env$comp
  lapply(1:max(comps), function(i) id[comps == i])
}

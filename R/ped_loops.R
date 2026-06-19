#' Pedigree loops
#'
#' Functions for identifying, breaking and restoring loops in pedigrees.
#'
#' In pedtools and other pedsuite packages, pedigree loops are usually handled
#' under the hood by the functions described here, without the need for explicit
#' action from end users. When a ped object `x` is created, an internal routine
#' detects if the pedigree contains loops, in which case `x$UNBROKEN_LOOPS` is
#' set to TRUE.
#'
#' In cases with complex inbreeding, it can be instructive to plot the pedigree
#' after breaking the loops. Duplicated individuals are plotted with suitable
#' labels (see examples).
#'
#' The function `findLoopBreakers()` searches the marriage graph for cycles.
#' For each cycle it chooses a suitable loop breaker and one child indicating
#' which full-sibship should be moved.
#'
#' The function `breakLoops()` then makes one founder copy of each selected
#' loop breaker. For each pair of loop breaker and child, all children with the
#' same parents as the indicated child are moved from the original parent to
#' the copy. Repeating this removes the cycles, while recording the original
#' and copied individuals in `x$LOOP_BREAKERS` so that `tieLoops()` can restore
#' the original pedigree. In the current implementation, only nonfounders can
#' act as loop breakers.
#'
#' Loop breakers may also be supplied by the user. A vector gives the loop
#' breakers themselves; for each one, a full-sibship is chosen automatically.
#' Repeated values are allowed and select successive sibships of the same
#' individual. Alternatively, `loopBreakers` may be a two-column matrix or data
#' frame whose first column gives loop breakers and whose second column gives
#' a child in the sibship to be moved.
#'
#' @param x A [ped()] object.
#' @param loopBreakers Either NULL, resulting in automatic selection of loop
#'   breakers, a vector indicating the individuals to be used as loop breakers,
#'   or a two-column matrix/data frame giving loop breakers and children.
#' @param verbose A logical indicating whether to print progress messages.
#' @param score Optional preference scores used by `findLoopBreakers()`. Higher
#'   values are preferred. If named, the names are interpreted as individual
#'   labels. Otherwise the vector must have length `pedsize(x)`.
#' @param errorIfFail A logical. If TRUE, an error is raised if loop breaking
#'   fails. If FALSE, the pedigree is returned unchanged.
#' @param allowFounder,allowRepeated Logicals, temporarily set to FALSE.
#'
#' @return For `breakLoops()`, a `ped` object in which the indicated loop
#'   breakers are duplicated. The returned object has a non-null
#'   `LOOP_BREAKERS` entry: a matrix with the internal IDs of the originals in
#'   the first column and the copies in the second. If loop breaking fails,
#'   then depending on `errorIfFail`, either an error is raised or the input
#'   pedigree is returned unchanged.
#'
#'   For `tieLoops()`, a `ped` object in which any duplicated individuals
#'   listed in `x$LOOP_BREAKERS` are merged back into their originals. For any
#'   ped object `x`, the call `tieLoops(breakLoops(x))` should return `x`.
#'
#'   For `findLoopBreakers()`, a two-column character matrix giving suitable
#'   loop breakers in the first column and split children in the second.
#'
#' @examples
#'
#' x = cousinPed(1, child = TRUE)
#' plot(breakLoops(x))
#'
#' # Pedigree with marriage loop: Double first cousins
#' y = doubleCousins(1, 1, child = TRUE)
#' y2 = breakLoops(y)
#' plot(y2)
#'
#' # Alternatively, first find loop breakers, then break loops with them
#' lb = findLoopBreakers(y)
#' y3 = breakLoops(y, lb)
#' stopifnot(identical(y2, y3))
#'
#' # Different set of loop breakers enforced
#' y4 = breakLoops(y, 6:7)
#' plot(y4)
#'
#' @importFrom utils tail
#' @importFrom stats ave
#' @export
breakLoops = function(x, loopBreakers = NULL, verbose = TRUE, errorIfFail = TRUE,
                      score = NULL, allowFounder = FALSE, allowRepeated = FALSE) {

  if(isFALSE(x$UNBROKEN_LOOPS) || is.singleton(x)) {
    if(verbose) {
      if(is.null(x$LOOP_BREAKERS)) message("No loops to break")
      else message("No further loops to break")
    }
    return(x)
  }

  labs = x$ID
  n = length(labs)
  auto = length(loopBreakers) == 0

  if(auto) {
    plan = findLoopBreakers(x, score = score, errorIfFail = errorIfFail,
                            allowFounder = allowFounder, allowRepeated = allowRepeated)
    if(!nrow(plan)) {
      if(errorIfFail) stop2("No loop breakers found")
      if(verbose) message("No loop breakers found - returning unchanged")
      return(x)
    }
  }
  else if(is.matrix(loopBreakers) || is.data.frame(loopBreakers)) {
    plan = as.matrix(loopBreakers)
  }
  else {
    lbInt = internalID(x, loopBreakers) # just to check validity
    if(any(lbInt %in% leaves(x, internal = TRUE)))
      stop2("Leaves cannot be loop breakers: ", .myintersect(leaves(x), loopBreakers))

    # Occurrence number of each lb among repeated values
    occ = ave(lbInt, lbInt, FUN = seq_along)

    chs = vapply(seq_along(lbInt), function(i) {
      sibs = children(x, lbInt[i], internal = TRUE, bySpouse = TRUE)
      if(occ[i] > length(sibs))
        stop2("Loop breaker repeated too many times: ", labs[lbInt[i]])
      utils::tail(sibs[[occ[i]]], 1)
    }, integer(1))

    plan = cbind(loopBreaker = labs[lbInt], child = labs[chs])
  }

  breakers = internalID(x, as.character(plan[, 1]))
  splitChild = internalID(x, as.character(plan[, 2]))

  oldpedm = as.matrix(x)
  attrs = attributes(oldpedm)
  oldLB = attrs$LOOP_BREAKERS
  oldOrigs  = oldLB[, 1] %||% integer(0)
  oldCopies = oldLB[, 2] %||% integer(0)

  if(any(breakers %in% oldCopies))
    stop2("Loop-breaker copies cannot be copied again")

  if(!allowFounder && any(breakers %in% founders(x, internal = TRUE)))
    stop2("Founder loop breakers require `allowFounder = TRUE`")

  if(!allowRepeated && anyDuplicated.default(c(oldOrigs, breakers)))
    stop2("Repeated loop breakers require `allowRepeated = TRUE`")

  if(verbose)
    cat("Loop breakers:", toString(labs[breakers]), "\n")

  # Copy selected rows
  cnt = tabulate(breakers, n)
  all_rows = rep.int(seq_len(n), 1L + cnt)
  new_rows = duplicated.default(all_rows)
  pedm = oldpedm[all_rows, , drop = FALSE]

  # New matrix IDs for the copies
  ord = order(breakers, seq_along(breakers))
  dups = integer(length(breakers))
  dups[ord] = n + seq_along(breakers)

  pedm[new_rows, 1] = n + seq_along(breakers)
  pedm[new_rows, 2:3] = 0L

  # Move selected child groups to their copies
  for(i in seq_along(breakers)) {
    r = match(splitChild[i], pedm[, 1])
    fa = pedm[r, 2]
    mo = pedm[r, 3]
    b = breakers[i]

    if(fa != b && mo != b)
      stop2("Individual ", labs[b], " is not a parent of ", labs[splitChild[i]])

    same = pedm[, 2] == fa & pedm[, 3] == mo
    if(fa == b)
      pedm[same, 2] = dups[i]
    if(mo == b)
      pedm[same, 3] = dups[i]
  }

  # Labels
  oldCount = tabulate(oldOrigs, n)
  copyNo = integer(length(breakers))
  copyNo[ord] = sequence.default(cnt[cnt > 0L])
  copyNo = copyNo + oldCount[breakers]

  newlabs = attrs$LABELS[all_rows]
  newlabs[match(dups, pedm[, 1])] = paste0(strrep("=", copyNo), labs[breakers])
  attrs$LABELS = newlabs

  # Loop breaker matrix
  newLB = rbind(oldLB, cbind(orig = breakers, copy = dups))
  newLB[] = match(newLB, pedm[, 1])
  attrs$LOOP_BREAKERS = newLB

  # Restore as ped object
  newx = restorePed(pedm, attrs = attrs)

  if(auto && isTRUE(newx$UNBROKEN_LOOPS)) {
    mess = "Loop breaking incomplete"
    if(errorIfFail) stop2(mess)
    if(verbose) message(mess)
  }

  newx
}

#' @export
#' @rdname breakLoops
tieLoops = function(x, verbose = TRUE) {
  dups = x$LOOP_BREAKERS
  if (is.null(dups) || nrow(dups) == 0) {
    if(verbose) cat("No loops to tie\n")
    return(x)
  }
  if (any(dups > pedsize(x)))
    stop2("Something's wrong - duplicated individual is out of range: ",
          dups[dups > pedsize(x)])

  origs = dups[, 1]
  copies = dups[, 2]

  oldpedm = as.matrix(x)
  attrs = attributes(oldpedm)

  # Remove copy labels
  attrs$LABELS = attrs$LABELS[-copies]
  attrs$LOOP_BREAKERS = NULL

  # Discard the duplicated rows
  newpedm = oldpedm[-copies, ]

  # Restore wrong parents
  wrong = match(newpedm[,2:3], copies, nomatch = 0)
  newpedm[,2:3][wrong > 0] = origs[wrong]

  restorePed(newpedm, attrs = attrs)
}


#' @export
#' @rdname breakLoops
findLoopBreakers = function(x, score = NULL, errorIfFail = TRUE,
                            allowFounder = FALSE, allowRepeated = FALSE) {

  empty = matrix(character(0), ncol = 2,
                 dimnames = list(NULL, c("lb", "child")))

  if(!isTRUE(x$UNBROKEN_LOOPS))
    return(empty)

  labs = x$ID
  n = length(labs)
  pref = numeric(n)

  if(!is.null(score)) {
    if(!is.null(names(score)))
      pref[internalID(x, names(score))] = score
    else if(length(score) == n)
      pref = score
    else
      stop2("`score` must have length pedsize(x), or be named")
  }

  oldLB = x$LOOP_BREAKERS
  oldOrigs = if(is.null(oldLB)) integer(0) else oldLB[, 1]
  oldCopies = if(is.null(oldLB)) integer(0) else oldLB[, 2]

  FOU = founders(x, internal = TRUE)

  idx = nonfounders(x, internal = TRUE)
  fidx = x$FIDX[idx]
  midx = x$MIDX[idx]
  breakers = splitChild = integer(length(idx))
  nb = 0L

  repeat {
    N = max(idx, fidx, midx)
    couples = fidx * (N + 1L) + midx
    g = .marriageGraphEdges(idx, fidx, midx, reduced = TRUE)

    loopRows = .findGraphCycle(cbind(g[, 1:2, drop = FALSE], seq_len(nrow(g))))
    if(!length(loopRows))
      break

    ed = g[as.integer(loopRows), , drop = FALSE]
    b = ed[, "label"]
    toNuc = match(ed[, "to"], couples, nomatch = 0L)

    cand = which(toNuc > 0L & b <= n & b %notin% oldCopies)
    s = pref[b[cand]]
    cand = cand[!is.na(s) & s > -Inf]

    usedBreakers = c(oldOrigs, breakers[seq_len(nb)])

    # TODO: These are intended to be temporary
    if(!allowFounder)
      cand = cand[b[cand] %notin% FOU]

    if(!allowRepeated)
      cand = cand[b[cand] %notin% usedBreakers]

    if(!length(cand)) {
      if(errorIfFail)
        stop2("The selected individuals cannot break all loops")
      return(empty)
    }

    used = b[cand] %in% usedBreakers
    pick = cand[order(-pref[b[cand]], !used)[1L]]
    j = toNuc[pick]

    nb = nb + 1L
    breakers[nb] = b[pick]
    splitChild[nb] = idx[j]

    # Move this child group to a temporary copy
    same = couples == couples[j]
    dup = n + nb
    if(fidx[j] == breakers[nb])
      fidx[same] = dup
    if(midx[j] == breakers[nb])
      midx[same] = dup
  }

  if(nb == 0L)
    return(empty)

  k = seq_len(nb)
  cbind(lb = labs[breakers[k]], child = labs[splitChild[k]])
}



# Internal method for creating marriage graph edges.
# NB: All internal/numerical labels.
.marriageGraphEdges = function(idx, fidx, midx, reduced = FALSE) {

  # Nuclear-family nodes, encoded safely above individual IDs
  N = max(c(idx, fidx, midx))
  couples = fidx*(N + 1L) + midx

  # Incoming edges: spouses to nuclear-family nodes
  keepNuc = !duplicated.default(couples)
  nuc = couples[keepNuc]
  fa = fidx[keepNuc]
  mo = midx[keepNuc]
  nonself = fa != mo
  spou = c(fa, mo[nonself])
  marIn = cbind(from = spou, to = c(nuc, nuc[nonself]), label = spou)

  # Outgoing edges: nuclear-family nodes to children
  marOut = cbind(from = couples, to = idx, label = idx)

  if(reduced) {
    # Leaves cannot contribute to cycles
    isNonleaf = idx %in% c(fidx, midx)
    marOut = marOut[isNonleaf, , drop = FALSE]

    # Contract individuals with only one spouse/family connection
    ismono = tabulate(spou, nbins = N) == 1L
    monos = which(ismono)
    monoOut = match(monos, marOut[, "to"], nomatch = 0L)
    hit = monoOut > 0L

    # Replace child monos with their marriage node
    if(any(hit)) {
      childMono = monos[hit]
      marOut[monoOut[hit], "to"] = marIn[match(childMono, marIn[, "from"]), "to"]
    }

    # Remove spouse edge
    marIn = marIn[!ismono[marIn[, "from"]], , drop = FALSE]
  }

  rbind(marIn, marOut)
}

# Finds a cycle in a graph given as an edge matrix.
.findGraphCycle = function(g, edgelabels = NULL) {
  # g: edge matrix with columns 'from', 'to', 'label' (optional)
  # Output: vector of vertices (or edges if 'label') forming a cycle, or NULL

  if(!length(g))
    return(NULL)

  .from = g[, 1L]
  .to = g[, 2L]
  edgelabels = edgelabels %||% if(ncol(g) > 2L) as.character(g[, 3L])
  hasLabel = !is.null(edgelabels)

  # Compact node labels to 1:nn; original labels are kept in nodes
  nodes = unique.default(c(.from, .to))
  .from = match(.from, nodes)
  .to = match(.to, nodes)
  nn = length(nodes)

  # Undirected edges
  a = .from
  b = .to
  sw = a > b
  a[sw] = .to[sw]
  b[sw] = .from[sw]
  edgeKey = a*(nn + 1L) + b

  # Parallel edges are 2-cycles in this reduced multigraph
  dup = duplicated.default(edgeKey)

  if(any(dup)) {
    j = which(dup)[1]
    i = match(edgeKey[j], edgeKey)
    return(if(hasLabel) edgelabels[c(i, j)] else as.character(nodes[c(a[j], b[j])]))
  }

  if(nn < 3L)
    return(NULL)

  # Adjacency list, undirected
  adjList = split.default(c(.to, .from), factor(c(.from, .to), levels = seq_len(nn)))

  env = new.env(parent = emptyenv())
  env$visited = logical(nn)
  env$prev = integer(nn)
  env$cycle = NULL

  # Depth first traversal, keeping track of path
  DFS = function(i, prev = 0L) {
    env$visited[i] = TRUE
    env$prev[i] = prev

    neigh = adjList[[i]]
    if(prev)
      neigh = neigh[neigh != prev]

    for(j in neigh) {
      if(env$visited[j]) {
        cyc = i
        x = i
        while(x != j) {
          x = env$prev[x]
          cyc = c(x, cyc)
        }
        env$cycle = cyc
        return(TRUE)
      }

      if(DFS(j, prev = i))
        return(TRUE)
    }

    FALSE
  }

  # Search all components
  for(v in seq_len(nn)) {
    if(!env$visited[v] && DFS(v))
      break
  }

  if(is.null(env$cycle))
    return(NULL)

  # Cycle given as compact vertex indices
  cycleV = env$cycle

  # If no edge label column, return original vertex labels
  if(!hasLabel)
    return(as.character(nodes[cycleV]))

  # Match cycle edges back to original edge labels
  a = cycleV
  b = c(cycleV[-1], cycleV[1])
  sw = a > b
  tmp = a[sw]
  a[sw] = b[sw]
  b[sw] = tmp

  edgelabels[match(a*(nn + 1L) + b, edgeKey)]
}


# Marriage graph of a pedigree. Not used, but useful for debugging and visualization.
# g = marriageGraph(x)
# gr = igraph::graph_from_edgelist(g[, 1:2], directed = FALSE)
# plot(gr) # further args needed to make it nice
marriageGraph = function(x, reduced = FALSE) {
  # Index of nonfounders
  idx = nonfounders(x, internal = TRUE)
  fidx = x$FIDX[idx]
  midx = x$MIDX[idx]

  g = .marriageGraphEdges(idx, fidx, midx, reduced = reduced)

  glab = as.character(g)
  dim(glab) = dim(g)

  N = max(c(idx, fidx, midx)) # NB: same as in .marriageGraphEdges()
  marnode = g > N
  fid = g[marnode] %/% (N + 1L)
  mid = g[marnode] %% (N + 1L)
  glab[marnode] = paste(fid, mid, sep = "+")

  glab
}



# LEGACY CODE ---------------------------------------------------------------------------------

#' Inbreeding loops
#'
#' Returns a list of all inbreeding loops in the pedigree. Note that this may be large
#' (and slow) in heavily inbred pedigrees. For example, `fullSibMating(10)` has 8144
#' loops. For finding loop breakers and breaking loops, see `findLoopBreakers()` and
#' `breakLoops()`, which are much more efficient for large pedigrees, and do not require
#' enumerating all loops.
#'
#' @param x A `ped` object.
#'
#' @return A list of inbreeding loops, where each loop is a list with entries
#' `top`, `bottom`, `pathA` and `pathB`. The `top` and `bottom` entries give
#' the IDs of the individuals at the top and bottom of the loop, while `pathA`
#' and `pathB` give the IDs of the individuals along the two paths from top to
#' bottom (excluding top and bottom themselves).
#'
#' @seealso [findLoopBreakers()], [breakLoops()]
#'
#' @examples
#' x = cousinPed(1, child = TRUE)
#' loops = inbreedingLoops(x)
#'
#' @export
inbreedingLoops = function(x) {

  if(!is.ped(x))
    stop2("`x` must be a ped object")

  n0 = pedsize(x)
  idmap = seq_len(n0)

  # NB: Very slow for large pedigrees. DFS probably more efficient.

  if(n0 > 30) {
    # Identify unneeded sib leaves
    pars = x$FIDX*(n0 + 1L) + x$MIDX
    sibships = unname(split.default(seq_along(pars), pars))
    sibships = sibships[lengths(sibships) > 1]

    lvsInt = leaves(x, TRUE)
    remov = unlist(lapply(sibships, function(s) {
      lvs = s %in% lvsInt
      if(all(lvs)) s[-1] else s[lvs]
    }))

    # Remove sib leaves and update idmap
    if(length(remov)) {
      keep = .mysetdiff(seq_len(n0), remov)
      x = x |> setMarkers(NULL) |> removeIndividuals(x$ID[remov], verbose = FALSE)
     idmap = idmap[keep]
    }
  }

  # Start main algorithm
  n = pedsize(x)
  dls = descentPaths(x, 1:n, internal = TRUE)
  dls = dls[lengths(dls) > 1]

  loops = lapply(dls, function(dl) {
    top = dl[[1]][1]
    pairs = .comb2(length(dl))
    loopsA = vector("list", length = nrow(pairs))
    i = 1
    for(p in 1:nrow(pairs)) {
      p1 = dl[[pairs[p, 1]]][-1]
      p2 = dl[[pairs[p, 2]]][-1]
      if (p1[1] == p2[1]) # skip if collapse (same member one step down)
        next
      inters = p1[match(p2, p1, 0L)] #intersecting path members, excluding id
      if(!length(inters))
        next
      bottom = inters[1]
      pathA = p1[seq_len(match(bottom, p1)-1)]  #without top and bottom. Seq_len to deal with the 1:0 problem.
      pathB = p2[seq_len(match(bottom, p2)-1)]
      loopsA[[i]] = list(top = top, bottom = bottom, pathA = pathA, pathB = pathB)
      i = i + 1
    }
    length(loopsA) = i - 1
    unique.default(loopsA)
  })

  loops = unlist(loops, recursive = FALSE)

  # Return original IDs
  lapply(loops, function(lo) {
    lo$top = idmap[lo$top]
    lo$bottom = idmap[lo$bottom]
    lo$pathA = idmap[lo$pathA]
    lo$pathB = idmap[lo$pathB]
    lo
  })
}



.breakLoops_old = function(x, loopBreakers = NULL, verbose = TRUE, errorIfFail = TRUE) {

  if (isFALSE(x$UNBROKEN_LOOPS) || is.singleton(x)) {
    if (verbose) {
      if(is.null(x$LOOP_BREAKERS)) message("No loops to break")
      else message("No further loops to break")
    }
    return(x)
  }

  auto = is.null(loopBreakers)
  if (auto) {
    loopBreakers = .findLoopBreakers_old(x)
    if (length(loopBreakers) == 0) {
      if (verbose)
        message("Marriage loops detected, trying different selection method")
      loopBreakers = .findLoopBreakers2_old(x, errorIfFail = errorIfFail)
    }
  }
  if (!length(loopBreakers)) {
    mess = "Loop breaking unsuccessful"
    if(errorIfFail) stop2(mess)
    else {
      if(verbose) message(mess, " - returning unchanged")
      return(x)
    }
  }

  LABS = labels(x)

  # Convert to internal IDs and sort (don't skip this)
  loopBreakers = .mysortInt(internalID(x, loopBreakers))

  FOU = founders(x, internal = TRUE)
  FOU_LB = intersect(loopBreakers, FOU)
  if (length(FOU_LB) > 0) {
    mess = "Breaking loops at founders is not implemented"
    if(errorIfFail) stop2(mess)
    else {
      if(verbose) message(mess, " - returning unchanged")
      return(x)
    }
  }

  if (verbose)
    cat("Loop breakers:", toString(LABS[loopBreakers]), "\n")

  ### Old ped data
  oldpedm = as.matrix(x)
  n = pedsize(x)

  ### New ped matrix
  # Setup for creating pedm by replicating lb rows
  all_rows = rep.int(1:n, times = 1 + (1:n %in% loopBreakers))
  new_rows = duplicated.default(all_rows)

  # Dummy numerical IDs of the new duplicated indivs.
  # NB: These are inserted in 'new_rows' positions, hence disrupts 1,2,3,...
  # Therefore they will always be changed in restorePed().
  dups = n + seq_along(loopBreakers)

  # Create it
  pedm = oldpedm[all_rows, ]
  pedm[new_rows, 1] = dups
  pedm[new_rows, 2:3] = 0

  # Change original loop breakers occurring in FIDX and MIDX
  wrong = match(pedm[, 2:3], loopBreakers, nomatch = 0)
  pedm[, 2:3][wrong > 0] = dups[wrong]

  ### Modify labels
  attrs = attributes(oldpedm)  #all attributes except 'dim'
  newlabs = attrs$LABELS[all_rows]
  newlabs[new_rows] = paste0("=", newlabs[new_rows])
  attrs$LABELS = newlabs

  ### Loop breaker matrix.
  # Previous loopBreakers must be fixed!
  if(!is.null(old_lb_matr <- attrs$LOOP_BREAKERS))
    old_lb_matr[] = match(old_lb_matr, pedm[, 1])

  # New rows. (Remember internal numbering, i.e. indices.)
  new_lb_matr = cbind(orig = which(new_rows) - 1, copy = which(new_rows))

  # Append the new lb matrix to the old
  attrs$LOOP_BREAKERS = rbind(old_lb_matr, new_lb_matr)

  ### Create new ped
  newx = restorePed(pedm, attrs = attrs)
  if (auto)
    return(.breakLoops_old(newx, verbose = verbose, errorIfFail = errorIfFail))
  newx
}



# NB: Legacy version of findLoopBreakers; only identifies and breaks inbreeding loops
.findLoopBreakers_old = function(x) {
  loopdata = inbreedingLoops(x)

  # write each loop as vector excluding top/bottom
  loops = lapply(loopdata, function(lo) c(lo$pathA, lo$pathB))
  bestbreakers = numeric()
  while (length(loops) > 0) {
    # add the individual occurring in most loops
    best = which.max(tabulate(unlist(loops)))
    bestbreakers = c(bestbreakers, best)
    loops = loops[sapply(loops, function(vec) best %notin% vec)]
  }
  labs = labels(x)
  labs[bestbreakers]
}

# Identifies and breaks marriage loops (was called after all inbreeding loops were broken)
.findLoopBreakers2_old = function(x, errorIfFail = TRUE) {

  breakers = numeric()
  N = pedsize(x)
  idx  = nonfounders(x, internal = TRUE)
  fidx = x$FIDX[idx]
  midx = x$MIDX[idx]
  nonf = as.character(idx)

  while (TRUE) {
    g = .marriageGraphEdges(idx, fidx, midx, reduced = TRUE)
    loop = .findGraphCycle(g)
    if(length(loop) == 0)
      break

    goodBreakers = .myintersect(loop, nonf)
    if(length(goodBreakers) == 0) {
      if(errorIfFail) stop("\
This pedigree requires founders as loop breakers, which is not implemented in pedtools yet.\
Please contact package maintainer if this is important to you.", call. = FALSE)
      else return()
    }

    b = as.numeric(goodBreakers[1])
    breakers = c(breakers, b)
    N = N+1
    fidx[fidx == b] = N
    midx[midx == b] = N
  }

  labs = labels(x)
  labs[breakers]
}

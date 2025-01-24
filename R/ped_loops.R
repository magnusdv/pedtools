#' Pedigree loops
#'
#' Functions for identifying, breaking and restoring loops in pedigrees.
#'
#' Pedigree loops are usually handled (by pedtools and related packages) under
#' the hood -- using the functions described here -- without the need for
#' explicit action from end users. When a ped object `x` is created, an internal
#' routine detects if the pedigree contains loops, in which case
#' `x$UNBROKEN_LOOPS` is set to TRUE.
#'
#' In cases with complex inbreeding, it can be instructive to plot the pedigree
#' after breaking the loops. Duplicated individuals are plotted with appropriate
#' labels (see examples).
#'
#' The function `breakLoops` breaks the loops of the input pedigree by
#' duplicating the *loop breakers*. These may be given by the user; otherwise
#' they are selected automatically. In the current implementation, only
#' nonfounders can act as loop breakers. For automatic selection of loop
#' breakers, `breakLoops` first calls `findLoopBreakers`, which identifies a set
#' of individuals breaking all *inbreeding loops* and breaks at the returned
#' individuals. If the resulting ped object still has loops, `findLoopBreakers2`
#' is called to handle *marriage loops*. In earlier versions of pedtools this
#' required the `igraph` package, but now uses a custom implementation using a
#' depth-first search algorithm to find a cycle in the marriage node graph of
#' the pedigree.
#'
#' @param x a [ped()] object.
#' @param loopBreakers either NULL (resulting in automatic selection of loop
#'   breakers) or a vector indicating the individuals to be used as loop
#'   breakers.
#' @param verbose a logical: Verbose output or not?
#' @param errorIfFail a logical: If TRUE an error is raised if the loop breaking
#'   is unsuccessful. If FALSE, the pedigree is returned unchanged.
#'
#' @return For `breakLoops`, a `ped` object in which the indicated loop breakers
#'   are duplicated. The returned object will also have a non-null
#'   `LOOP_BREAKERS` entry, namely a matrix with the IDs of the original loop
#'   breakers in the first column and the duplicates in the second. If loop
#'   breaking fails, then depending on `errorIfFail` either an error is raised,
#'   or the input pedigree is returned, still containing unbroken loops.
#'
#'   For `tieLoops`, a `ped` object in which any duplicated individuals (as
#'   given in the `x$LOOP_BREAKERS` entry) are merged. For any ped object `x`,
#'   the call `tieLoops(breakLoops(x))` should return `x`.
#'
#'   For `inbreedingLoops`, a list containing all inbreeding loops (not marriage
#'   loops) found in the pedigree. Each loop is represented as a list with
#'   elements `top`, `bottom`, `pathA` (individuals forming a path from top to
#'   bottom) and `pathB` (creating a different path from top to bottom, with no
#'   individuals in common with `pathA`). Note that the number of loops reported
#'   here counts all closed paths in the pedigree and will in general be larger
#'   than the genus of the underlying graph.
#'
#'   For `findLoopBreakers` and `findLoopBreakers2`, a vector of individual
#'   labels.
#'
#' @examples
#'
#' x = cousinPed(1, child = TRUE)
#' plot(breakLoops(x))
#'
#' # Pedigree with marriage loop: Double first cousins
#' y = doubleCousins(1, 1, child = TRUE)
#' findLoopBreakers(y) # --> 9
#' findLoopBreakers2(y) # --> 5 and 9
#' y2 = breakLoops(y)
#' plot(y2)
#'
#' # Or loop breakers chosen by user
#' y3 = breakLoops(y, 6:7)
#' plot(y3)
#'
#' @export
inbreedingLoops = function(x) { # CHANGE: pedigreeLoops changed name to inbreedingLoops
  n = pedsize(x)
  dls = descentPaths(x, 1:n, internal = TRUE)
  dls = dls[lengths(dls) > 1]

  loops = list()
  for (dl in dls) {
    top = dl[[1]][1]
    pairs = .comb2(length(dl))
    for (p in 1:nrow(pairs)) {
      p1 = dl[[pairs[p, 1]]][-1]
      p2 = dl[[pairs[p, 2]]][-1]
      if (p1[1] == p2[1]) # skip if collapse (same member one step down)
        next
      inters = p1[match(p2, p1, 0L)] #intersecting path members, excluding id
      if (!length(inters))
        next
      bottom = inters[1]
      pathA = p1[seq_len(match(bottom, p1)-1)]  #without top and bottom. Seq_len to deal with the 1:0 problem.
      pathB = p2[seq_len(match(bottom, p2)-1)]
      loops = c(loops, list(list(top = top, bottom = bottom, pathA = pathA, pathB = pathB)))
    }
  }
  unique(loops)
}

#' @export
#' @rdname inbreedingLoops
breakLoops = function(x, loopBreakers = NULL, verbose = TRUE, errorIfFail = TRUE) {

  if (isFALSE(x$UNBROKEN_LOOPS) || is.singleton(x)) {
    if (verbose) {
      if(is.null(x$LOOP_BREAKERS)) message("No loops to break")
      else message("No further loops to break")
    }
    return(x)
  }

  auto = is.null(loopBreakers)
  if (auto) {
    loopBreakers = findLoopBreakers(x)
    if (length(loopBreakers) == 0) {
      if (verbose)
        message("Marriage loops detected, trying different selection method")
      loopBreakers = findLoopBreakers2(x, errorIfFail = errorIfFail)
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

  # Change original loop breakers occuring in FIDX and MIDX
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
    return(breakLoops(newx, verbose = verbose, errorIfFail = errorIfFail))
  newx
}

#' @export
#' @rdname inbreedingLoops
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
#' @rdname inbreedingLoops
findLoopBreakers = function(x) {
  loopdata = inbreedingLoops(x)
  # write each loop as vector excluding top/bottom
  loops = lapply(loopdata, function(lo) c(lo$pathA, lo$pathB))
  bestbreakers = numeric()
  while (length(loops) > 0) {
    # add the individual occurring in most loops
    best = which.max(tabulate(unlist(loops)))
    bestbreakers = c(bestbreakers, best)
    loops = loops[sapply(loops, function(vec) !best %in% vec)]
  }
  labs = labels(x)
  labs[bestbreakers]
}

#' @export
#' @rdname inbreedingLoops
findLoopBreakers2 = function(x, errorIfFail = TRUE) {

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


### Jan 2025: Methods for finding and breaking marriage loops
marriageGraph = function(x, reduced = FALSE) {
  # Index of nonfounders
  idx = nonfounders(x, internal = TRUE)
  fidx = x$FIDX[idx]
  midx = x$MIDX[idx]

  g = .marriageGraphEdges(idx, fidx, midx, reduced = reduced)

  glab = as.character(g)
  dim(glab) = dim(g)
  marnodes = g[g > 1000]
  fid = marnodes %/% 1000
  mid = marnodes %% 1000
  glab[g > 1000] = paste(fid, mid, sep = "+")
  glab
}

# Internal method for creating marriage graph edges.
# NB: All internal/numerical labels.
.marriageGraphEdges = function(idx, fidx, midx, reduced = FALSE) {
  couples = 1000*fidx + midx

  # Incoming edges: From each spouse to their marriage node
  dups = duplicated.default(couples)
  spou = c(fidx[!dups], midx[!dups])
  marIn = cbind(from = spou, to = rep(couples[!dups], 2), label = spou)

  # Outgoing edges, from marriage nodes to children
  marOut = cbind(from = couples, to = idx, label = idx)

  if(reduced) {
    # Skip edges to leaves
    isnonleave = idx %in% c(fidx,midx)
    marOut = marOut[isnonleave, , drop = FALSE]

    # Skip ID nodes with only 1 spouse
    ismono = tabulate(spou) == 1
    monos = which(ismono)

    monoidxOut = match(monos, marOut[, 2], nomatch = 0)
    childMono = monos[monoidxOut > 0]

    # Replace child monos with their marriage node
    marOut[monoidxOut[monoidxOut > 0], 2] = marIn[match(childMono, marIn[, 1]), 2]

    # Remove spouse edge
    marIn = marIn[!ismono[marIn[,1]], , drop = FALSE]
  }

  rbind(marIn, marOut)
}

# Finds a cycle in a graph given as an edge matrix.
.findGraphCycle = function(g) {
  # g: edge matrix with columns 'from', 'to', 'label' (optional)
  # Output: vector of vertices (or edges if 'label') forming a cycle, or NULL

  mode(g) = "character"
  .from = g[,1]
  .to = g[,2]
  nodes = unique.default(c(.from, .to))

  # Adjacency list (undirected!)
  adjList = lapply(nodes, function(i) c(.from[.to == i], .to[.from == i]))
  names(adjList) = nodes

  env = new.env()
  env$visited = logical(length(adjList)) |> .setnames(nodes)
  env$prev  = rep(NA, length(adjList)) |> .setnames(nodes)
  env$cycle = NULL

  # Depth first traversal, keeping track of path
  DFS = function(i, prev = "") {
    env$visited[i] = TRUE
    env$prev[i]  = prev
    neigh = adjList[[i]] |> .mysetdiff(prev)
    for(j in neigh) {
      if(env$visited[j]) {
        # Reconstruct cycle
        cyc = x = i
        while(x != j)
          cyc = c(x <- env$prev[x], cyc)
        env$cycle = as.character(cyc)
        return(TRUE)
      }
      else if(DFS(j, prev = i))
        return(TRUE)
    }
    FALSE
  }

  # Start DFS from first node
  DFS(names(adjList)[length(nodes)])

  if(is.null(env$cycle))
    return(NULL)

  # Cycle given as vector of vertices
  cycleV = env$cycle

  # If no edge label column, return vertices
  if(ncol(g) == 2)
    return(cycleV)

  # Matrix of edge labels
  nn = length(nodes)
  edgeLabels = matrix(NA_character_, nrow = nn, ncol = nn, dimnames = list(nodes, nodes))
  edgeLabels[g[, 1:2, drop = FALSE]] = g[,3]
  edgeLabels[g[, 2:1, drop = FALSE]] = g[,3]

  # Return edge labels of cycle
  unique.default(edgeLabels[cbind(cycleV, c(cycleV[-1], cycleV[1]))])
}

#' Pedigree loops
#'
#' Functions for identifying, breaking and restoring loops in pedigrees.
#'
#' Pedigree loops are usually handled (by pedtools and related packages) under
#' the hood - using the functions described here - without need for explicit
#' action from end users. When a ped object `x` is created, an internal routine
#' detects if the pedigree contains loops, in which case `x$UNBROKEN_LOOPS` is
#' set to TRUE.
#'
#' In cases with complex inbreeding, it can be instructive to plot the pedigree
#' after breaking the loops. Duplicated individuals are plotted with appropriate
#' labels (see examples).
#'
#' The function `findLoopBreakers` identifies a set of individuals breaking all
#' inbreeding loops, but not marriage loops. These require more machinery for
#' efficient detection, and pedtools does this is a separate function,
#' `findLoopBreakers2`, utilizing methods from the `igraph` package. Since this
#' is rarely needed for most users, `igraph` is not imported when loading
#' pedtools, only when `findLoopBreakers2` is called.
#'
#' In practice, `breakLoops` first calls `findLoopBreakers` and breaks at the
#' returned individuals. If the resulting ped object still has loops,
#' `findLoopBreakers2` is called to break any marriage loops.
#'
#' @param x a [ped()] object.
#' @param loop_breakers either NULL (resulting in automatic selection of loop
#'   breakers) or a numeric containing IDs of individuals to be used as loop
#'   breakers.
#' @param verbose a logical: Verbose output or not?
#' @param errorIfFail a logical: If TRUE an error is raised if the loop breaking
#'   is unsuccessful. If FALSE, the pedigree is returned unchanged.
#' @return For `breakLoops`, a `ped` object in which the indicated loop breakers
#'   are duplicated. The returned object will also have a non-null
#'   `loop_breakers` entry, namely a matrix with the IDs of the original loop
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
#'   For `findLoopBreakers` and `findLoopBreakers2`, a numeric vector of
#'   individual ID's.
#' @author Magnus Dehli Vigeland
#'
#' @examples
#'
#' x = cousinPed(1, child = TRUE)
#' plot(breakLoops(x))
#'
#' # Pedigree with marriage loop: Double first cousins
#' if(requireNamespace("igraph", quietly = TRUE)) {
#'   y = doubleCousins(1, 1, child = TRUE)
#'   findLoopBreakers(y) # --> 5
#'   findLoopBreakers2(y) # --> 5 and 6
#'   y2 = breakLoops(y)
#'   plot(y2)
#'
#'   # Or loop breakers chosen by user
#'   y3 = breakLoops(y, c(3, 7))
#'   plot(y3)
#' }
#'
#' @export
inbreedingLoops = function(x) { # CHANGE: pedigreeLoops changed name to inbreedingLoops
  n = pedsize(x)
  dls = .descentPaths(x, 1:n, internal=TRUE)
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
breakLoops = function(x, loop_breakers = NULL, verbose = TRUE, errorIfFail = TRUE) {
  if (isFALSE(x$UNBROKEN_LOOPS) || is.singleton(x)) {
    if (verbose) {
      if(is.null(x$LOOP_BREAKERS)) message("No loops to break")
      else message("No further loops to break")
    }
    return(x)
  }

  auto = is.null(loop_breakers)
  if (auto) {
    loop_breakers = findLoopBreakers(x)
    if (length(loop_breakers) == 0) {
      if (verbose)
        message("Marriage loops detected, trying different selection method")
      loop_breakers = findLoopBreakers2(x, errorIfFail = errorIfFail)
    }
  }
  if (!length(loop_breakers)) {
    mess = "Loop breaking unsuccessful"
    if(errorIfFail) stop2(mess)
    else {
      if(verbose) message(mess, " - returning unchanged")
      return(x)
    }
  }

  LABS = labels(x)

  # Convert to internal IDs and sort (don't skip this)
  loop_breakers = .mysortInt(internalID(x, loop_breakers))

  FOU = founders(x, internal=T)
  FOU_LB = intersect(loop_breakers, FOU)
  if (length(FOU_LB) > 0) {
    mess = "Breaking loops at founders is not implemented"
    if(errorIfFail) stop2(mess)
    else {
      if(verbose) message(mess, " - returning unchanged")
      return(x)
    }
  }

  if (verbose)
    cat("Loop breakers:", toString(LABS[loop_breakers]), "\n")

  ### Old ped data
  oldpedm = as.matrix(x)  #data.frame(x, famid=T, missing=0)
  n = pedsize(x)

  ### New ped matrix
  # Setup for creating pedm by replicating lb rows
  all_rows = rep.int(1:n, times = 1 + (1:n %in% loop_breakers))
  new_rows = duplicated.default(all_rows)

  # Dummy numerical IDs of the new duplicated indivs.
  # NB: These are inserted in 'new_rows' positions, hence disrupts 1,2,3,...
  # Therefore they will always be changed in restorePed().
  dups = n + seq_along(loop_breakers)

  # Create it
  pedm = oldpedm[all_rows, ]
  pedm[new_rows, 1] = dups
  pedm[new_rows, 2:3] = 0

  # Change original loop breakers occuring in FIDX and MIDX
  wrong = match(pedm[,2:3], loop_breakers, nomatch=0)
  pedm[,2:3][wrong > 0] = dups[wrong]

  ### Modify labels
  attrs = attributes(oldpedm)  #all attributes except 'dim'
  newlabs = attrs$LABELS[all_rows]
  newlabs[new_rows] = paste0("=", newlabs[new_rows])
  attrs$LABELS = newlabs

  ### Loop breaker matrix.
  # Previous loop_breakers must be fixed!
  if(!is.null(old_lb_matr <- attrs$LOOP_BREAKERS))
    old_lb_matr[] = match(old_lb_matr, pedm[,1])

  # New rows. (Remember internal numbering, i.e. indices.)
  new_lb_matr = cbind(orig=which(new_rows)-1, copy=which(new_rows))

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
tieLoops = function(x, verbose=TRUE) {
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
  wrong = match(newpedm[,2:3], copies, nomatch=0)
  newpedm[,2:3][wrong > 0] = origs[wrong]

  restorePed(newpedm, attrs = attrs)
}

#' @export
#' @rdname inbreedingLoops
findLoopBreakers = function(x) {
  loopdata = inbreedingLoops(x)
  # write each loop as vector exluding top/bottom
  loops = lapply(loopdata, function(lo) c(lo$pathA, lo$pathB))
  bestbreakers = numeric()
  while (length(loops) > 0) {
    # add the individual occuring in most loops
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
  if (!requireNamespace("igraph", quietly = TRUE)) {
    warning("This pedigree has marriage loops. The package 'igraph' must be installed for automatic selection of loop breakers.\n")
    return(numeric())
  }

  breakers = numeric()
  N = pedsize(x)

  ped2edge = function(id, fid, mid) {
    # input: ped-kolonner UTEN founder-rader
    couples = paste(fid, mid, sep = "+")
    dups = duplicated.default(couples)
    edge.children = cbind(couples, id)
    edge.marriage_F = cbind(fid, couples)[!dups, ]
    edge.marriage_M = cbind(mid, couples)[!dups, ]
    rbind(edge.marriage_F, edge.marriage_M, edge.children)
  }

  NONFOU = nonfounders(x, internal=T)
  id = NONFOU
  fid = x$FIDX[NONFOU]
  mid = x$MIDX[NONFOU]
  nonf = as.character(NONFOU)

  while (T) {
    g = igraph::graph_from_edgelist(ped2edge(id, fid, mid))
    loop = igraph::girth(g)$circle
    if (length(loop) == 0)
      break
    good.breakers = intersect(loop$name, nonf)

    if (length(good.breakers) == 0) {
      if(errorIfFail) stop("\
This pedigree requires founders as loop breakers, which is not implemented in pedtools yet.\
Please contact magnusdv at medisin.uio.no if this is important to you.", call.=FALSE)
      else return()
    }

    b = as.numeric(good.breakers[1])
    breakers = c(breakers, b)
    N = N+1
    fid[fid == b] = N
    mid[mid == b] = N
  }
  labs = labels(x)
  labs[breakers]
}


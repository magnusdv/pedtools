#' Pedigree subgroups
#'
#' A collection of utility functions for identifying pedigree members with
#' certain properties.
#'
#' @param x A [ped()] object or a list of such.
#' @param id,ids A character (or coercible to such) with one or several ID
#'   labels.
#' @param inclusive A logical indicating whether an individual should be counted
#'   among his or her own ancestors/descendants
#' @param internal A logical indicating whether `id` (or `ids`) refers to the
#'   internal order.
#' @param degree,removal Non-negative integers.
#' @param half a logical or NA. If TRUE (resp. FALSE), only half (resp. full)
#'   siblings/cousins/nephews/nieces are returned. If NA, both categories are
#'   included.
#'
#' @return The functions `ancestors(x, id)` and `descendants(x, id)` return a
#'   vector containing the IDs of all ancestors (resp. descendants) of the
#'   individual `id` within the pedigree `x`. If `inclusive = TRUE`, `id` is
#'   included in the output.
#'
#'   For `commonAncestors(x, ids)` and `commonDescendants(x, ids)`, a vector
#'   containing the IDs of common ancestors to all of `ids`.
#'
#'   The functions `founders`, `nonfounders`, `males`, `females`, `leaves` each
#'   return a vector containing the IDs of all pedigree members with the wanted
#'   property. (Recall that a founder is a member without parents in the
#'   pedigree, and that a leaf is a member without children in the pedigree.)
#'
#'   The functions `father`, `mother`, `cousins`, `grandparents`,
#'   `nephews_nieces`, `children`, `parents`, `siblings`, `spouses`,
#'   `unrelated`, each returns a vector containing the IDs of all pedigree
#'   members having the specified relationship with `id`.
#'
#' @author Magnus Dehli Vigeland
#'
#' @examples
#'
#' x = ped(id = 2:9,
#'          fid = c(0,0,2,0,4,4,0,2),
#'          mid = c(0,0,3,0,5,5,0,8),
#'          sex = c(1,2,1,2,1,2,2,2))
#'
#' spouses(x, id = 2) # 3, 8
#' children(x, 2)     # 4, 9
#' descendants(x, 2)  # 4, 6, 7, 9
#' siblings(x, 4)     # 9 (full or half)
#' unrelated(x, 4)    # 5, 8
#' father(x, 4)       # 2
#' mother(x, 4)       # 3
#'
#' siblings(x, 4, half = FALSE) # none
#' siblings(x, 4, half = TRUE)  # 9
#'
#' leaves(x)          # 6, 7, 9
#' founders(x)        # 2, 3, 5, 8
#'
#' @name ped_subgroups
NULL


#' @rdname ped_subgroups
#' @export
founders = function(x, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    return(unname(unlist(lapply(x, founders, internal = FALSE))))
  }

  isFOU = x$FIDX == 0
  if (internal) which(isFOU) else labels.ped(x)[isFOU]
}

#' @rdname ped_subgroups
#' @export
nonfounders = function(x, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    return(unname(unlist(lapply(x, nonfounders, internal = FALSE))))
  }

  isNF = x$FIDX > 0
  if(internal) which(isNF) else labels.ped(x)[isNF]
}

#' @rdname ped_subgroups
#' @export
leaves = function(x, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    return(unname(unlist(lapply(x, leaves, internal = FALSE))))
  }

  lvs = if(is.singleton(x)) 1L else (1:pedsize(x))[-c(x$FIDX, x$MIDX)]
  if(internal) lvs else labels.ped(x)[lvs]
}

#' @rdname ped_subgroups
#' @export
males = function(x, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    return(unname(unlist(lapply(x, males, internal = FALSE))))
  }

  m = x$SEX == 1
  if(internal) which(m) else labels.ped(x)[m]
}

#' @rdname ped_subgroups
#' @export
females = function(x, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    return(unname(unlist(lapply(x, females, internal = FALSE))))
  }

  f = x$SEX == 2
  if(internal) which(f) else labels.ped(x)[f]
}

#' @rdname ped_subgroups
#' @export
typedMembers = function(x, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    return(unname(unlist(lapply(x, typedMembers, internal = FALSE))))
  }

  nMark = nMarkers(x)
  labs = x$ID
  if(nMark == 0)
    return(if(internal) integer(0) else character(0))

  allelematrix = unlist(x$MARKERS)
  typed = .rowSums(allelematrix, m = length(labs), n = 2*nMark) > 0

  # dim(allelematrix) = c(pedsize(x), 2*nMark)
  # typed = rowSums(allelematrix) > 0

  if(internal) which(typed) else labs[typed]
}


#' @rdname ped_subgroups
#' @export
untypedMembers = function(x, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    return(unname(unlist(lapply(x, untypedMembers, internal = FALSE))))
  }

  nMark = nMarkers(x)
  labs = x$ID
  if(nMark == 0)
    return(if(internal) seq_along(labs) else labs)

  allelematrix = unlist(x$MARKERS)
  untyped = .rowSums(allelematrix, m = length(labs), n = 2*nMark) == 0

  # dim(allelematrix) = c(pedsize(x), 2*nMark)
  # untyped = rowSums(allelematrix) == 0

  if(internal) which(untyped) else labs[untyped]
}

#' @rdname ped_subgroups
#' @export
father = function(x, id, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    comp = getComponent(x, id, checkUnique = TRUE, errorIfUnknown = TRUE)
    return(father(x[[comp]], id, internal = FALSE))
  }

  if(!internal)
    id = internalID(x, id)

  fa = x$FIDX[id]
  if(internal) fa else labels.ped(x)[fa]
}

#' @rdname ped_subgroups
#' @export
mother = function(x, id, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    comp = getComponent(x, id, checkUnique = TRUE, errorIfUnknown = TRUE)
    return(mother(x[[comp]], id, internal = FALSE))
  }

  if(!internal)
    id = internalID(x, id)

  mo = x$MIDX[id]
  if(internal) mo else labels.ped(x)[mo]
}

#' @rdname ped_subgroups
#' @export
children = function(x, id, internal = FALSE) {
  if(length(id) != 1)
    stop2("`id` must have length 1")

  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    comp = getComponent(x, id, checkUnique = TRUE, errorIfUnknown = TRUE)
    return(children(x[[comp]], id, internal = FALSE))
  }

  if(!internal)
    id = internalID(x, id)

  ch = (x$FIDX == id | x$MIDX == id)
  if(internal) which(ch) else labels.ped(x)[ch]
}

#' @rdname ped_subgroups
#' @export
offspring = children

#' @rdname ped_subgroups
#' @export
spouses = function(x, id, internal = FALSE) {
  if(length(id) != 1)
    stop2("`id` must have length 1")

  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    comp = getComponent(x, id, checkUnique = TRUE, errorIfUnknown = TRUE)
    return(spouses(x[[comp]], id, internal = FALSE))
  }

  # Returns a vector containing all individuals sharing offspring with <id>.
  if(!internal)
    id = internalID(x, id)

  spous = switch(x$SEX[id] + 1,
                c(x$MIDX[x$FIDX == id], x$FIDX[x$MIDX == id]), # sex = 0
                x$MIDX[x$FIDX == id],                        # sex = 1
                x$FIDX[x$MIDX == id])                        # sex = 2
  spous_uniq = unique.default(spous)
  if(internal) spous_uniq else labels.ped(x)[spous_uniq]
}


#' @rdname ped_subgroups
#' @export
unrelated = function(x, id, internal = FALSE) {
  if(length(id) != 1)
    stop2("`id` must have length 1")

  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    comp = getComponent(x, id, checkUnique = TRUE, errorIfUnknown = TRUE)
    unr = unrelated(x[[comp]], id, internal = FALSE)

    # Add indivs from all other comps
    unr = c(unr, unname(unlist(labels(x[-comp]))))
    return(unr)
  }

  if(!internal)
    id = internalID(x, id)

  ancs = ancestors(x, id, inclusive = TRUE, internal = TRUE)
  rel = lapply(ancs, function(a) descendants(x, a, inclusive = TRUE, internal = TRUE))
  unrel = setdiff(1:pedsize(x), unlist(rel))

  if(internal) unrel else labels.ped(x)[unrel]
}


#' @rdname ped_subgroups
#' @export
parents = function(x, id, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    comp = getComponent(x, id, checkUnique = TRUE, errorIfUnknown = TRUE)
    return(parents(x[[comp]], id, internal = FALSE))
  }

  if(!internal)
    id = internalID(x, id)

  par = c(x$FIDX[id], x$MIDX[id])
  if(internal) par else labels.ped(x)[par]
}

#' @rdname ped_subgroups
#' @export
grandparents = function(x, id, degree = 2, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    comp = getComponent(x, id, checkUnique = TRUE, errorIfUnknown = TRUE)
    return(grandparents(x[[comp]], id, degree = degree, internal = FALSE))
  }

  if(!internal)
    id = internalID(x, id)

  nextgen = id
  for(i in seq_len(degree))
    nextgen = c(x$FIDX[nextgen], x$MIDX[nextgen])

  if(internal) nextgen else labels.ped(x)[nextgen]
}

#' @rdname ped_subgroups
#' @export
siblings = function(x, id, half = NA, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")
    comp = getComponent(x, id, checkUnique = TRUE, errorIfUnknown = TRUE)
    return(siblings(x[[comp]], id, half = half, internal = FALSE))
  }

  if(!internal)
    id = internalID(x, id)

  fa = x$FIDX[id]
  mo = x$MIDX[id]
  if(fa == 0 && mo == 0)
    return(if(internal) integer(0) else character(0))

  samefather = x$FIDX == fa
  samemother = x$MIDX == mo

  sib =
    if(isTRUE(half)) xor(samefather, samemother)   # half only
    else if(isFALSE(half)) samefather & samemother # full only
    else if(is.na(half)) samefather | samemother   # either
  sib[id] = FALSE
  if(internal) which(sib) else labels.ped(x)[sib]
}


# TODO: Review this before re-export
cousins = function(x, id, degree = 1, removal = 0, half = NA, internal = FALSE) {
  if (!internal)  id = internalID(x, id)
  gp = grandparents(x, id, degree = degree, internal = TRUE)
  gp = gp[gp > 0]
  if(length(gp) == 0)
    return(if(internal) integer(0) else character(0))

  uncles = unique.default(unlist(lapply(gp, function(a)
    siblings(x, a, half = half, internal = TRUE))))

  cous = uncles
  for (i in seq_len(degree + removal))
    cous = unique.default(unlist(lapply(cous, children, x = x, internal = TRUE)))

  if (internal) cous else labels.ped(x)[cous]
}

#' @rdname ped_subgroups
#' @export
nephews_nieces = function(x, id, removal = 1, half = NA, internal = FALSE) {
    cousins(x, id, degree = 0, removal = removal, half = half, internal = internal)
}

#' @rdname ped_subgroups
#' @export
ancestors = function(x, id, inclusive = FALSE, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")

    comps = getComponent(x, id, checkUnique = TRUE, errorIfUnknown = TRUE)
    ancList = lapply(unique.default(comps), function(co) {
      idsComp = id[comps == co]
      ancestors(x[[co]], idsComp, inclusive = inclusive, internal = FALSE)
    })
    return(unlist(ancList))
  }

  # climbs upwards storing parents iteratively. (Not documented: Accepts id of length > 1)
  if(!internal)
    id = internalID(x, id)

  FIDX = x$FIDX
  MIDX = x$MIDX
  ancest = if(inclusive) id else integer(0)

  up1 = c(FIDX[id], MIDX[id])
  up1 = up1[up1 > 0]
  while (length(up1)) {
    ancest = c(ancest, up1)
    up1 = c(FIDX[up1], MIDX[up1])
    up1 = up1[up1 > 0]
  }
  ancest = .mysortInt(unique.default(ancest))
  if(internal) ancest else labels.ped(x)[ancest]
}

#' @rdname ped_subgroups
#' @export
commonAncestors = function(x, ids, inclusive = FALSE, internal = FALSE) {
  if(length(ids) < 2)
    stop2("Argument `ids` must have length at least 2")

  anc = ancestors(x, ids[1], inclusive = inclusive, internal = internal)
  for(id in ids[-1]) {
    if(length(anc) == 0)
      break
    newanc = ancestors(x, id, inclusive = inclusive, internal = internal)
    anc = .myintersect(anc, newanc)
  }

  anc
}

#' @rdname ped_subgroups
#' @export
descendants = function(x, id, inclusive = FALSE, internal = FALSE) {
  if(is.pedList(x)) {
    if(internal)
      stop2("Argument `internal` cannot be TRUE when `x` is a pedlist")

    comps = getComponent(x, id, checkUnique = TRUE, errorIfUnknown = TRUE)
    ancList = lapply(unique.default(comps), function(co) {
      idsComp = id[comps == co]
      descendants(x[[co]], idsComp, inclusive = inclusive, internal = FALSE)
    })
    return(unlist(ancList))
  }

  if(!internal)
    id = internalID(x, id)

  FIDX = x$FIDX
  MIDX = x$MIDX
  desc = if(inclusive) id else integer()

  nextoffs = id
  while(length(nextoffs)) {
      nextoffs = which(FIDX %in% nextoffs | MIDX %in% nextoffs)
      desc = c(desc, nextoffs)
  }
  desc = .mysortInt(unique.default(desc))
  if(internal) desc else labels.ped(x)[desc]
}

#' @rdname ped_subgroups
#' @export
commonDescendants = function(x, ids, inclusive = FALSE, internal = FALSE) {
  if(length(ids) < 2)
    stop2("Argument `ids` must have length at least 2")

  desc = descendants(x, ids[1], inclusive = inclusive, internal = internal)
  for(id in ids[-1]) {
    if(length(desc) == 0)
      break
    newdesc = descendants(x, id, inclusive = inclusive, internal = internal)
    desc = .myintersect(desc, newdesc)
  }

  desc
}


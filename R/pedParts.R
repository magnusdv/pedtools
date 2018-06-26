#' Pedigree subsets
#'
#' Utility functions for 'ped' objects, mainly for extracting various
#' pedigree information.
#'
#' @param x A [ped()] object.
#' @param id A single ID label (numeric or character).
#' @param internal A logical indicating whether 'id' refers to the internal order.
#' @param degree,removal Non-negative integers.
#' @param half a logical or NA. If TRUE (resp FALSE), only half (resp. full)
#' siblings/cousins/nephews/nieces are returned. If NA, both categories are
#' inclucded.
#'
#' @return For `ancestors(x,id)`, a vector containing the ID's of all
#' ancestors of the individual `id`.  For `descendants(x,id)`, a
#' vector containing the ID's of all descendants (i.e. children, grandchildren,
#' a.s.o.) of individual `id`.
#'
#' The functions `cousins`, `grandparents`, `nephews_nieces`,
#' `children`, `parents`, `siblings`, `spouses`,
#' `unrelated`, each returns an integer vector containing the ID's of all
#' pedigree members having the specified relationship with `id`.
#'
#' For `leaves`, a vector of the ID labels of pedigree members without
#' children.
#'
#' @author Magnus Dehli Vigeland
#'
#' @examples
#'
#' x = ped(id=2:9, fid=c(0,0,2,0,4,4,0,2), mid=c(0,0,3,0,5,5,0,8),
#'         sex=c(1,2,1,2,1,2,2,2))
#' stopifnot(setequal(spouses(x, 2), c(3,8)),
#'           setequal(children(x, 2), c(4,9)),
#'           setequal(descendants(x, 2), c(4,6,7,9)),
#'           setequal(leaves(x), c(6,7,9)))
#'
#' @name pedParts
NULL


#' @rdname pedParts
#' @export
founders = function(x, internal = FALSE) {
  if (internal) x$FOUNDERS else x$LABELS[x$FOUNDERS]
}

#' @rdname pedParts
#' @export
nonfounders = function(x, internal = FALSE) {
  if (internal) x$NONFOUNDERS else x$LABELS[x$NONFOUNDERS]
}

#' @rdname pedParts
#' @export
children = function(x, id, internal = FALSE) {
    if (!internal) id = internalID(x, id)
    offs_int = (x$FID == id | x$MID == id)

    if (internal) which(offs_int) else x$LABELS[offs_int]
}

#' @rdname pedParts
#' @export
offspring = children

#' @rdname pedParts
#' @export
spouses = function(x, id, internal = FALSE) {
  # Returns a vector containing all individuals sharing offspring with <id>.
  if (!internal)  id = internalID(x, id)
  spous = switch(x$SEX[id] + 1,
                c(x$MID[x$FID == id], x$FID[x$MID == id]), # sex = 0
                x$MID[x$FID == id],                        # sex = 1
                x$FID[x$MID == id])                        # sex = 2
  spous_uniq = unique.default(spous)
  if (internal) spous_uniq else x$LABELS[spous_uniq]
}


#' @rdname pedParts
#' @export
unrelated = function(x, id, internal = FALSE) {
  if (!internal)  id = internalID(x, id)
  ancs = c(id, ancestors(x, id))
    rel = unique.default(unlist(lapply(ancs, function(a) c(a, descendants(x, a, internal = FALSE)))))
    unrel = setdiff(x$LABELS, rel)
    if (internal) internalID(x, unrel) else unrel
}


#' @rdname pedParts
#' @export
leaves = function(x, internal = FALSE) {
  leaves_int = setdiff(x$ID, c(x$FID, x$MID))
  if (internal) leaves_int else x$LABELS[leaves_int]
}

#' @rdname pedParts
#' @export
parents = function(x, id, internal = FALSE) {
  if (!internal) id = internalID(x, id)
  parents_int = c(x$FID[id], x$MID[id])
  if (internal) parents_int else x$LABELS[parents_int]
}

#' @rdname pedParts
#' @export
grandparents = function(x, id, degree = 2, internal = FALSE) {
  if (!internal)  id = internalID(x, id)

  nextgen = id
  for (i in seq_len(degree)) nextgen = c(x$FID[nextgen], x$MID[nextgen])
  if (internal) nextgen else x$LABELS[nextgen]
}

#' @rdname pedParts
#' @export
siblings = function(x, id, half = NA, internal = FALSE) {
  if (!internal)  id = internalID(x, id)
  if (id %in% x$FOUNDERS) return(numeric(0))
  fa = x$FID[id]
  mo = x$MID[id]

  samefather = x$FID == fa
  samemother = x$MID == mo
  sib_int =
    if (isTRUE(half)) samefather | samemother
    else if (isFALSE(half)) xor(samefather, samemother)
    else if(is.na(half)) samefather & samemother
  sib_int[id] = FALSE
  if (internal) which(sib_int) else x$LABELS[sib_int]
}

#' @rdname pedParts
#' @export
cousins = function(x, id, degree = 1, removal = 0, half = NA, internal = FALSE) {
  if (!internal)  id = internalID(x, id)
  gp = grandparents(x, id, degree = degree, internal = TRUE)
  uncles = unique.default(unlist(lapply(gp, siblings, x = x, half = half, internal = TRUE)))
  cous = uncles
  for (i in seq_len(degree + removal))
    cous = unique.default(unlist(lapply(cous, children, x = x, internal = TRUE)))

  if (internal) cous else x$LABELS[cous]
}

#' @rdname pedParts
#' @export
nephews_nieces = function(x, id, removal = 1, half = NA, internal = FALSE) {
    cousins(x, id, degree = 0, removal = removal, half = half, internal = internal)
}

#' @rdname pedParts
#' @export
ancestors = function(x, id, internal = FALSE) {
  # climbs upwards storing parents iteratively. (Not documented: Accepts id of length > 1)
  if (!internal)  id = internalID(x, id)
  FID = x$FID
  MID = x$MID
  ancest = numeric(0)
  up1 = c(FID[id], MID[id])
  up1 = up1[up1 > 0]
  while (length(up1)) {
    ancest = c(ancest, up1)
    up1 = c(FID[up1], MID[up1])
    up1 = up1[up1>0]
  }
  ancest = sort.int(unique.default(ancest))
  if (internal) ancest else x$LABELS[ancest]
}


#' @rdname pedParts
#' @export
descendants = function(x, id, internal = FALSE) {
  if (!internal)  id = internalID(x, id)

  FID = x$FID
  MID = x$MID

  desc = numeric()
  nextoffs = id
  while(length(nextoffs)) {
      nextoffs = which(FID %in% nextoffs | MID %in% nextoffs)
      desc = c(desc, nextoffs)
  }
  desc = sort.int(unique.default(desc))
  if (internal) desc else x$LABELS[desc]
}


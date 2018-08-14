#' Pedigree subsets
#'
#' Utility functions for 'ped' objects, mainly for extracting various pedigree
#' information.
#'
#' @param x A [ped()] object.
#' @param id A single ID label (coercible to character).
#' @param internal A logical indicating whether 'id' refers to the internal
#'   order.
#' @param degree,removal Non-negative integers.
#' @param half a logical or NA. If TRUE (resp FALSE), only half (resp. full)
#'   siblings/cousins/nephews/nieces are returned. If NA, both categories are
#'   inclucded.
#'
#' @return For `ancestors(x,id)`, a vector containing the IDs of all ancestors
#'   of the individual `id`.  For `descendants(x,id)`, a vector containing the
#'   IDs of all descendants (i.e. children, grandchildren, a.s.o.) of
#'   individual `id`.
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
#' x = ped(id=2:9, fid=c(0,0,2,0,4,4,0,2), mid=c(0,0,3,0,5,5,0,8),
#'         sex=c(1,2,1,2,1,2,2,2))
#' stopifnot(setequal(spouses(x, 2), c(3,8)),
#'           setequal(children(x, 2), c(4,9)),
#'           setequal(descendants(x, 2), c(4,6,7,9)),
#'           setequal(leaves(x), c(6,7,9)))
#'
#' @name ped_subsets
NULL


#' @rdname ped_subsets
#' @export
founders = function(x, internal = FALSE) {
  is_fou = x$FIDX == 0
  if (internal) which(is_fou) else labels(x)[is_fou]
}

#' @rdname ped_subsets
#' @export
nonfounders = function(x, internal = FALSE) {
  is_nonfou = x$FIDX > 0
  if (internal) which(is_nonfou) else labels(x)[is_nonfou]
}

#' @rdname ped_subsets
#' @export
leaves = function(x, internal = FALSE) {
  leaves_int = (1:pedsize(x))[-c(x$FIDX, x$MIDX)]
  if (internal) leaves_int else labels(x)[leaves_int]
}

#' @rdname ped_subsets
#' @export
males = function(x, internal = FALSE) {
  m = x$SEX == 1
  if (internal) which(m) else labels(x)[m]
}

#' @rdname ped_subsets
#' @export
females = function(x, internal = FALSE) {
  f = x$SEX == 2
  if (internal) which(f) else labels(x)[f]
}

#' @rdname ped_subsets
#' @export
typedMembers = function(x, internal = FALSE) {
  if (nMarkers(x) == 0)
    return(if(internal) numeric(0) else character(0))
  allelematrix = do.call(cbind, x$markerdata)
  emptyrows = rowSums(allelematrix != 0) == 0
  if(internal) which(!emptyrows) else labels(x)[!emptyrows]
}

#' @rdname ped_subsets
#' @export
untypedMembers = function(x, internal = FALSE) {
  if (nMarkers(x) == 0)
    return(if(internal) seq_len(pedsize(x)) else labels(x))
  allelematrix = do.call(cbind, x$markerdata)
  emptyrows = rowSums(allelematrix != 0) == 0
  if(internal) which(emptyrows) else labels(x)[emptyrows]
}

#' @rdname ped_subsets
#' @export
father = function(x, id, internal = FALSE) {
  if (!internal) id = internalID(x, id)
  fa = x$FIDX[id]
  if (internal) fa else labels(x)[fa]
}

#' @rdname ped_subsets
#' @export
mother = function(x, id, internal = FALSE) {
  if (!internal) id = internalID(x, id)
  mo = x$MIDX[id]
  if (internal) mo else labels(x)[mo]
}

#' @rdname ped_subsets
#' @export
children = function(x, id, internal = FALSE) {
    if (!internal) id = internalID(x, id)
    offs_int = (x$FIDX == id | x$MIDX == id)

    if (internal) which(offs_int) else labels(x)[offs_int]
}

#' @rdname ped_subsets
#' @export
offspring = children

#' @rdname ped_subsets
#' @export
spouses = function(x, id, internal = FALSE) {
  # Returns a vector containing all individuals sharing offspring with <id>.
  if (!internal)  id = internalID(x, id)
  spous = switch(x$SEX[id] + 1,
                c(x$MIDX[x$FIDX == id], x$FIDX[x$MIDX == id]), # sex = 0
                x$MIDX[x$FIDX == id],                        # sex = 1
                x$FIDX[x$MIDX == id])                        # sex = 2
  spous_uniq = unique.default(spous)
  if (internal) spous_uniq else labels(x)[spous_uniq]
}


#' @rdname ped_subsets
#' @export
unrelated = function(x, id, internal = FALSE) {
  if (!internal)  id = internalID(x, id)
  ancs = c(id, ancestors(x, id))
    rel = unique.default(unlist(lapply(ancs, function(a) c(a, descendants(x, a, internal = FALSE)))))
    unrel = setdiff(labels(x), rel)
    if (internal) internalID(x, unrel) else unrel
}


#' @rdname ped_subsets
#' @export
parents = function(x, id, internal = FALSE) {
  if (!internal) id = internalID(x, id)
  parents_int = c(x$FIDX[id], x$MIDX[id])
  if (internal) parents_int else labels(x)[parents_int]
}

#' @rdname ped_subsets
#' @export
grandparents = function(x, id, degree = 2, internal = FALSE) {
  if (!internal)  id = internalID(x, id)

  nextgen = id
  for (i in seq_len(degree)) nextgen = c(x$FIDX[nextgen], x$MIDX[nextgen])
  if (internal) nextgen else labels(x)[nextgen]
}

#' @rdname ped_subsets
#' @export
siblings = function(x, id, half = NA, internal = FALSE) {
  if (!internal)  id = internalID(x, id)
  fa = x$FIDX[id]
  mo = x$MIDX[id]
  if (fa==0 && mo==0) return(numeric(0))

  samefather = x$FIDX == fa
  samemother = x$MIDX == mo
  sib_int =
    if (isTRUE(half)) samefather | samemother
    else if (isFALSE(half)) xor(samefather, samemother)
    else if(is.na(half)) samefather & samemother
  sib_int[id] = FALSE
  if (internal) which(sib_int) else labels(x)[sib_int]
}

#' @rdname ped_subsets
#' @export
cousins = function(x, id, degree = 1, removal = 0, half = NA, internal = FALSE) {
  if (!internal)  id = internalID(x, id)
  gp = grandparents(x, id, degree = degree, internal = TRUE)
  uncles = unique.default(unlist(lapply(gp, siblings, x = x, half = half, internal = TRUE)))
  cous = uncles
  for (i in seq_len(degree + removal))
    cous = unique.default(unlist(lapply(cous, children, x = x, internal = TRUE)))

  if (internal) cous else labels(x)[cous]
}

#' @rdname ped_subsets
#' @export
nephews_nieces = function(x, id, removal = 1, half = NA, internal = FALSE) {
    cousins(x, id, degree = 0, removal = removal, half = half, internal = internal)
}

#' @rdname ped_subsets
#' @export
ancestors = function(x, id, internal = FALSE) {
  # climbs upwards storing parents iteratively. (Not documented: Accepts id of length > 1)
  if (!internal)  id = internalID(x, id)
  FIDX = x$FIDX
  MIDX = x$MIDX
  ancest = numeric(0)
  up1 = c(FIDX[id], MIDX[id])
  up1 = up1[up1 > 0]
  while (length(up1)) {
    ancest = c(ancest, up1)
    up1 = c(FIDX[up1], MIDX[up1])
    up1 = up1[up1>0]
  }
  ancest = sort.int(unique.default(ancest))
  if (internal) ancest else labels(x)[ancest]
}


#' @rdname ped_subsets
#' @export
descendants = function(x, id, internal = FALSE) {
  if (!internal)  id = internalID(x, id)

  FIDX = x$FIDX
  MIDX = x$MIDX

  desc = numeric()
  nextoffs = id
  while(length(nextoffs)) {
      nextoffs = which(FIDX %in% nextoffs | MIDX %in% nextoffs)
      desc = c(desc, nextoffs)
  }
  desc = sort.int(unique.default(desc))
  if (internal) desc else labels(x)[desc]
}


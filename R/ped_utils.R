#' Pedigree utilities
#'
#' Various utility functions for `ped` objects
#'
#' The functions `hasUnbrokenLoops()`, `hasInbredFounders()` and
#' `hasSelfing()` allow as input either a single `ped` object or a list of
#' such. In the latter case each function returns TRUE if it is TRUE for any of
#' the components.
#'
#' @param x A `ped` object, or (in some functions - see Details) a list of such.
#' @param chromType Either "autosomal" (default) or "x".
#'
#' @return
#'
#' * `pedsize(x)` returns the number of pedigree members in `x`
#'
#' * `hasUnbrokenLoops(x)` returns TRUE if `x` has loops, otherwise FALSE. (No
#' computation is done here; the function simply returns the value of
#' `x$UNBROKEN_LOOPS`).
#'
#' * `hasInbredFounders(x)` returns TRUE is founder inbreeding is specified
#' for `x` and at least one founder has positive inbreeding coefficient. See
#' [founderInbreeding()] for details.
#'
#' * `hasSelfing(x)` returns TRUE if the pedigree contains selfing events. This
#' is recognised by father and mother begin equal for some child. (Note that for
#' this to be allowed, the gender code of the parent must be 0.)
#'
#' * `hasCommonAncestor(x)` computes a logical matrix `A` whose entry `A[i,j]`
#' is TRUE if pedigree members i and j have a common ancestor in `x`, and FALSE
#' otherwise. By convention, `A[i,i]` is TRUE for all i.
#'
#' * `subnucs(x)` returns a list of all nuclear sub-pedigrees of `x`, wrapped as
#' `nucleus` objects. Each nucleus is a list with entries `father`, `mother` and
#' `children`.
#'
#' * `peelingOrder(x)` calls `subnucs(x)` and extends each entry with a `link`
#' individual, indicating a member linking the nucleus to the remaining
#' pedigree. One application of this function is the fact that if _fails_ to
#' find a complete peeling order if and only if the pedigree has loops. (In fact
#' it is called each time a new `ped` object is created by [ped()] in order to
#' detect loops.) The main purpose of the function, however, is to prepare for
#' probability calculations in other packages, as e.g. in
#' `pedprobr::likelihood`.
#'
#' @examples
#' x = fullSibMating(1)
#' stopifnot(pedsize(x) == 6)
#' stopifnot(hasUnbrokenLoops(x))
#'
#' # All members have common ancestors except the grandparents
#' CA = hasCommonAncestor(x)
#' stopifnot(!CA[1,2], !CA[2,1], sum(CA) == length(CA) - 2)
#'
#' # Effect of breaking the loop
#' y = breakLoops(x)
#' stopifnot(!hasUnbrokenLoops(y))
#' stopifnot(pedsize(y) == 7)
#'
#' # A pedigree with selfing (note the necessary `sex = 0`)
#' z1 = singleton(1, sex = 0)
#' z2 = addChildren(z1, father = 1, mother = 1, nch = 1)
#' stopifnot(!hasSelfing(z1), hasSelfing(z2))
#'
#' # Nucleus sub-pedigrees
#' stopifnot(length(subnucs(z1)) == 0)
#' peelingOrder(cousinPed(1))
#'
#' @name ped_utils
NULL


#' @rdname ped_utils
#' @export
pedsize = function(x) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  length(x$ID)
}


#' @rdname ped_utils
#' @export
hasUnbrokenLoops = function(x) {
  if(is.pedList(x))
    return(any(vapply(x, hasUnbrokenLoops, logical(1))))

  isTRUE(x$UNBROKEN_LOOPS)
}


#' @rdname ped_utils
#' @export
hasInbredFounders = function(x, chromType = "autosomal") {
  if(is.pedList(x))
    return(any(vapply(x, hasInbredFounders, logical(1), chromType = chromType)))

  if(is.null(x$FOUNDER_INBREEDING))
    return(FALSE)

  finb = founderInbreeding(x, named = T, chromType = chromType)

  # If X: only females interesting (males are always 1)
  if(chromType == "x")
    finb = finb[getSex(x, names(finb)) == 2]

  any(finb > 0)
}


#' @rdname ped_utils
#' @export
hasSelfing = function(x) {
  if(is.pedList(x))
    return(any(vapply(x, hasSelfing, logical(1))))

  any(x$FIDX != 0 & x$FIDX == x$MIDX)
}


#' @rdname ped_utils
#' @export
hasCommonAncestor = function(x) {
  if(!is.ped(x))
    stop2("Input is not a `ped` object: ", x)

  n = pedsize(x)
  labs = labels(x)

  A = matrix(FALSE, ncol = n, nrow = n, dimnames = list(labs, labs))

  FOU = founders(x, internal = TRUE)
  for(i in FOU) {
    # vector of all descendants of i, including i
    desc = c(i, descendants(x, i, internal = TRUE))
    A[fast.grid(rep(list(desc), 2))] = TRUE
  }
  A
}

validate_sex = function(sex, nInd, zero_allowed = TRUE) {
  if(length(sex) == 0)
    stop2(sprintf("`%s` cannot be empty", deparse(substitute(sex))))
  if(length(sex) > nInd)
    stop2(sprintf("`%s` is longer than the number of individuals",
                  deparse(substitute(sex))))

  codes = if(zero_allowed) 0:2 else 1:2
  sex_int = suppressWarnings(as.integer(sex))
  ok = (sex == sex_int) & sex_int %in% codes
  if(!all(ok))
    stop2("Illegal gender code: ", unique(sex[!ok]),
          ".  [0 = NA; 1 = male; 2 = female]")
  invisible(rep_len(sex_int, nInd))
}



#' @rdname ped_utils
#' @export
subnucs = function(x) {
  if (is.singleton(x))
    return(list())

  n = pedsize(x)
  seqn = seq_len(n)

  FIDX = x[["FIDX"]]
  MIDX = x[["MIDX"]]

  # Indices of unique parent couples
  p_pairs_idx = seqn[FIDX + MIDX > 0 & !duplicated.default(FIDX*(n+1) + MIDX)]

  # List all nucs: Format = c(father, mother, children)
  lapply(rev(p_pairs_idx), function(j) {
    nuc = list(father = FIDX[j], mother = MIDX[j],
               children = seqn[FIDX == FIDX[j] & MIDX == MIDX[j]])
    class(nuc) = "nucleus"
    attr(nuc, "labels") = labels(x)
    nuc
  })
}


#' @rdname ped_utils
#' @export
peelingOrder = function(x) {
  # output: list of nuclear subfamilies. Format for each nuc:
  # list(father,mother,children,link), where link = 0 for the last nuc.
  nucs = subnucs(x)
  if(length(nucs) == 0)
    return(nucs)

  peeling = vector("list", Nnucs <- length(nucs))
  i = k = 1

  while (length(nucs)) {
    # Start searching for a nuc with only one link to others
    nuc = nucs[[i]]

    # Identify links to other remaining nucs
    if(length(nucs) > 1) {
      nucmembers = c(nuc$father, nuc$mother, nuc$children)
      links = nucmembers[nucmembers %in% unlist(nucs[-i], use.names = F)]
    }
    else
      links = 0 # if nuc is the last

    # If only one link: move nuc to peeling list, and proceed
    if (length(links) == 1) {
      nuc$link = links
      peeling[[k]] = nuc
      nucs[i] = NULL
      i = 1
      k = k+1
    }
    else {
      if (i == length(nucs)) { # LOOP! Include remaining nucs without 'link', and break.
        peeling[k:Nnucs] = nucs
        break
      }
      # Otherwise try next remaining nuc
      i = i+1
    }
  }
  peeling
}


#' S3 methods
#'
#' @param x An object
#' @param ... Not used
#' @export
print.nucleus = function(x, ...) {
  labs = attr(x, 'labels')
  fa = x$father
  mo = x$mother
  ch = x$children

  link = x$link
  if (is.null(link))
    linktext = "NO LINK ASSIGNED"
  else if (link == 0)
    linktext = "0 (end of chain)"
  else {
    who = if (link == fa) "father"
    else if (link == mo) "mother"
    else if (link %in% ch) "child"
    else stop2("Erroneous link individual (not a member of nucleus): ", link)
    linktext = sprintf("%s (%s)", labs[link], who)
  }
  cat(sprintf("Nucleus: Fa/mo/ch = %s/%s/%s\nLink individual: %s\n",
              labs[fa], labs[mo], paste(labs[ch], collapse = ","), linktext))
}


has_numlabs = function(x) {
  # Returns TRUE if the labels of x are coercible to integers
  labs = labels(x)
  numlabs = suppressWarnings(as.character(as.integer(labs)))
  isTRUE(all(labs == numlabs))
}


.generations = function(x) {
  FOU = founders(x, internal = T)
  max(lengths(unlist(.descentPaths(x, FOU, internal = T), recursive = F)))
}

# Utility function for generating numbered "NN" labels.
# Returns "NN_i" where i increments largest j occuring as NN_j, NN.j or NN-j in input.
nextNN = function(labs) { # labs a character vector
  NNs = grepl("^NN", labs)
  if(!any(NNs))
    return("NN_1")
  NNnum = suppressWarnings(as.numeric(sub("^NN[._-]?", "", labs[NNs])))
  if(all(is.na(NNnum)))
    return("NN_1")
  nextNNnum = max(NNnum, na.rm = T) + 1
  return(sprintf("NN_%d", nextNNnum))
}


.descentPaths = function(x, ids, internal = FALSE) {
  if (!internal) ids = internalID(x, ids)

  offs = lapply(1:pedsize(x), children, x = x, internal = TRUE)
  lapply(ids, function(id) {
    res = list(id)
    while (TRUE) {
      newoffs = offs[vapply(res, function(path) path[length(path)], 1)]
      if (length(unlist(newoffs)) == 0)
        break
      nextstep = lapply(1:length(res), function(r)
        if (length(newoffs[[r]]) == 0) res[r]
        else lapply(newoffs[[r]], function(kid) c(res[[r]], kid)))
      res = unlist(nextstep, recursive = FALSE)
    }
    if (!internal) {
      labs = labels(x)
      res = lapply(res, function(int_ids) labs[int_ids])
    }
    res
  })
}

#KAN KASTES:
.descentPaths2 = function(x, ids, internal = FALSE) {
  if (!internal) ids = internalID(x, ids)

  offs = lapply(1:pedsize(x), function(y) children(x, y, internal = TRUE))
  lapply(ids, function(id) {
    res = list(id)
    while (TRUE) {
      newoffs = offs[vapply(res, function(path) path[length(path)], 1)]
      if (length(unlist(newoffs)) == 0)
        break
      nextstep = lapply(1:length(res), function(r)
        if (length(newoffs[[r]]) == 0) res[r]
        else lapply(newoffs[[r]], function(kid) c(res[[r]], kid)))
      res = unlist(nextstep, recursive = FALSE)
    }
    if (!internal) {
      labs = labels(x)
      res = lapply(res, function(int_ids) labs[int_ids])
    }
    res
  })
}



#' Pedigree component
#'
#' Given a list of `ped` objects (called pedigree components), and a vector of
#' ID labels, find the index of the component holding each individual.
#' @param x A list of `ped` objects
#' @param ids A vector of ID labels (coercible to character)
#' @param checkUnique If TRUE an error is raised if any element of `ids` occurs
#'   more than once in `x`.
#'
#' @return An integer vector of the same length as `ids`, with NA entries where
#'   the corresponding label was not found in any of the components.
#'
#' @examples
#' x = list(nuclearPed(1), singleton(id = "A"))
#' getComponent(x, c(3, "A")) # = c(1, 2)
#'
#' @export
getComponent = function(x, ids, checkUnique = FALSE) {
  if(is.ped(x))
    x = list(x)
  else if(!is.pedList(x))
    stop2("Input is not a (list of) ped objects")

  # List labels of each component
  labList = lapply(x, labels)

  # A single vector with all labels
  labVec = unlist(labList)

  # Check for duplicates if indicated
  if(checkUnique) {
    v = labVec[labVec %in% ids]
    if(dup <- anyDuplicated.default(v))
      stop2("ID label is not unique: ", v[dup])
  }

  # Vector of same length as labVec, with component index for each member
  compi = rep(seq_along(labList), times = lengths(labList))

  # Match input ids against label vector
  idx = match(ids, labVec)

  # Return comp idx of the input ids, including NA if not present
  compi[idx]
}

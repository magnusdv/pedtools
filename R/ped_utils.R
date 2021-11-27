#' Pedigree utilities
#'
#' Various utility functions for `ped` objects.
#'
#' @param x A `ped` object, or (in some functions - see Details) a list of such.
#' @param maxOnly A logical, by default TRUE. (See Value.)
#' @param maxComp A logical, by default TRUE. (See Value.)
#' @param chromType Either "autosomal" (default) or "x".
#'
#' @return
#'
#' * `pedsize(x)` returns the number of pedigree members in each component of
#' `x`.
#'
#' * `generations(x)` by default returns the number of generations in `x`,
#' defined as the number of individuals in the longest line of parent-child
#' links. (Note that this definition is valid also if `x` has loops.) If
#' `maxOnly = FALSE`, the output is a named integer vector, showing the
#' generation number of each pedigree member. If `x` has multiple components,
#' the output depends on the parameter `maxComp`. If this is FALSE, the output
#' is a vector containing the result for each component. If TRUE (default), only
#' the highest number is returned.
#'
#' * `hasUnbrokenLoops(x)` returns TRUE if `x` has loops, otherwise FALSE. (No
#' computation is done here; the function simply returns the value of
#' `x$UNBROKEN_LOOPS`).
#'
#' * `hasInbredFounders(x)` returns TRUE is founder inbreeding is specified for
#' `x` and at least one founder has positive inbreeding coefficient. See
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
#' stopifnot(generations(x) == 3)
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
  if(is.ped(x))
    return(length(x$ID))
  else if (is.pedList(x))
    return(vapply(x, function(comp) length(comp$ID), 1L))
  stop2("Input to `pedsize()` must be a `ped` object or a list of such")
}

#' @rdname ped_utils
#' @export
generations = function(x, maxOnly = TRUE, maxComp = TRUE) {

  if(is.pedList(x)) {
    gens = sapply(x, generations)
    return(if(maxComp) max(gens) else gens)
  }

  xorig = x
  x = parentsBeforeChildren(x)

  FIDX = x$FIDX
  MIDX = x$MIDX
  NONFOU = nonfounders(x, internal = TRUE)
  N = pedsize(x)

  # Vector of (maximal) generation number of each ID: dp[i] = 1 + max(dp[parents])
  # NB: Requires "parentsBeforeChildren".
  dp = rep(1L, N)
  for(i in NONFOU)
    dp[i] = 1L + max(dp[c(FIDX[i], MIDX[i])])

  if(maxOnly)
    return(max(dp))

  names(dp) = x$ID
  dp[xorig$ID]
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

  finb = founderInbreeding(x, named = TRUE, chromType = chromType)

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
    stop2("Illegal sex: ", unique(sex[!ok]),
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
      links = nucmembers[nucmembers %in% unlist(nucs[-i], use.names = FALSE)]
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


hasNumLabs = function(x) {
  # Returns TRUE if the labels of x are coercible to integers
  labs = labels(x)
  numlabs = suppressWarnings(as.character(as.integer(labs)))
  isTRUE(all(labs == numlabs))
}


.generations = function(x) {
  FOU = founders(x, internal = TRUE)
  max(lengths(unlist(descentPaths(x, FOU, internal = TRUE), recursive = FALSE)))
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
  nextNNnum = max(NNnum, na.rm = TRUE) + 1
  return(sprintf("NN_%d", nextNNnum))
}


#' Pedigree component
#'
#' Given a list of `ped` objects (called pedigree components), and a vector of
#' ID labels, find the index of the component holding each individual.
#'
#' @param x A list of `ped` objects
#' @param ids A vector of ID labels (coercible to character)
#' @param checkUnique If TRUE an error is raised if any element of `ids` occurs
#'   more than once in `x`. Default: FALSE.
#' @param errorIfUnknown If TRUE, the function stops with an error if not all
#'   elements of `ids` are recognised as names of members in `x`. Default:
#'   FALSE.
#'
#' @return An integer vector of the same length as `ids`, with NA entries where
#'   the corresponding label was not found in any of the components.
#'
#' @seealso [internalID()]
#'
#' @examples
#' x = list(nuclearPed(1), singleton(id = "A"))
#' getComponent(x, c(3, "A"))
#'
#' @export
getComponent = function(x, ids, checkUnique = FALSE, errorIfUnknown = FALSE) {
  if(is.ped(x))
    x = list(x)
  else if(!is.pedList(x))
    stop2("Input is not a (list of) ped objects")

  # List labels of each component
  labList = labels(x)

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

  if(anyNA(idx) && errorIfUnknown)
    stop2("Unknown ID label: ", ids[is.na(idx)])

  # Return comp idx of the input ids, including NA if not present
  compi[idx]
}



topolOrder = function(x) {
  N = pedsize(x)

  env = new.env()
  logvec = logical(N)
  names(logvec) = x$ID

  isNonFou = logvec
  isNonFou[nonfounders(x)] = TRUE

  env$visited = logvec
  env$ord = character(N)
  env$k = 1

  dfs = function(x, id, env) {
    env$visited[id] = TRUE

    if(isNonFou[id]) {
      fa = father(x, id)
      if(!env$visited[fa])
        dfs(x, fa, env)
      mo = mother(x, id)
      if(!env$visited[mo])
        dfs(x, mo, env)
    }

    k = env$k
    env$ord[k] = id
    #message("Added ", id, " in position ", N-k)

    env$k = k + 1
  }

  lvs = sort.int(leaves(x), method = "quick")
  for(id in lvs)
    dfs(x, id, env)

  # Return ordered labels
  env$ord
}

reorderSimple = function(x, neworder, ids = NULL) {
  if(!is.ped(x))
    stop2("`reorderQuick` only works for `ped` objects")
  if(length(neworder) != length(x$ID))
    stop2("Wrong length of `neworder`")

  if(is.character(neworder))
    neworder = match(neworder, x$ID)

  ID = x$ID[neworder]
  FIDX = x$FIDX[neworder]
  MIDX = x$MIDX[neworder]
  SEX = x$SEX[neworder]

  nonf = FIDX > 0
  FIDX[nonf] = match(FIDX[nonf], neworder)
  MIDX[nonf] = match(MIDX[nonf], neworder)

  # Mask all but `ids` indivs
  if(!is.null(ids)) {
    mask = !ID %in% ids
    ID[mask] = "*"
    SEX[mask] = "*"
  }

  list(ID, FIDX, MIDX, SEX)
}


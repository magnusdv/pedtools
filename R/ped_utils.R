#' Pedigree utilities
#'
#' Various utility functions for `ped` objects
#'
#' @param x A `ped` object
#'
#' @return
#'
#' * `pedsize(x)` returns the number of pedigree members in `x`
#'
#' * `has_unbroken_loops(x)` returns TRUE if `x` has loops, otherwise FALSE. (No
#' computation is done here; the function simply returns the value of
#' `x$UNBROKEN_LOOPS`)
#'
#' * `has_inbred_founders(x)` returns TRUE is founder inbreeding is specified
#' for `x` AND at least one founder has positive inbreeding coefficient. See
#' [founder_inbreeding()] for details.
#'
#' * `has_selfing(x)` returns TRUE if the pedigree contains selfing events. This
#' is recognised by father and mother begin equal for some child. (Note that for
#' this to be allowed, the gender code of the parent must be 0.)
#'
#' * `has_common_ancestor(x)` computes a logical matrix `A` whose entry
#' `A[i,j]` is TRUE if pedigree members i and j have a common ancestor in `x`,
#' and FALSE otherwise. By convention, `A[i,i]` is TRUE for all i.
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
#' [pedprobr::likelihood()].
#'
#' @examples
#' x = fullSibMating(2)
#' stopifnot(pedsize(x) == 6)
#' stopifnot(has_unbroken_loops(x))
#'
#' # All members have common ancestors except the grandparents
#' CA = has_common_ancestor(x)
#' stopifnot(!CA[1,2], !CA[2,1], sum(CA) == length(CA) - 2)
#'
#' # Effect of breaking the loop
#' y = breakLoops(x)
#' stopifnot(!has_unbroken_loops(y))
#' stopifnot(pedsize(y) == 7)
#'
#' # A pedigree with selfing (note the neccessary `sex = 0`)
#' z1 = singleton(1, sex = 0)
#' z2 = addChildren(z1, father = 1, mother = 1, nch = 1)
#' stopifnot(!has_selfing(z1), has_selfing(z2))
#'
#' # Nucleus sub-pedigrees
#' stopifnot(length(subnucs(z1)) == 0)
#' peelingOrder(cousinsPed(1))
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
has_unbroken_loops = function(x) {
  isTRUE(x$UNBROKEN_LOOPS)
}


#' @rdname ped_utils
#' @export
has_inbred_founders = function(x) {
  finb = x$FOUNDER_INBREEDING
  !is.null(finb) && any(finb > 0)
}


#' @rdname ped_utils
#' @export
has_selfing = function(x) {
  any(x$FIDX != 0 & x$FIDX == x$MIDX)
}


#' @rdname ped_utils
#' @export
has_common_ancestor = function(x) {
  n = pedsize(x)
  labs = labels(x)

  A = matrix(FALSE, ncol=n, nrow=n, dimnames=list(labs, labs))

  FOU = founders(x, internal=TRUE)
  for(i in FOU) {
    # vector of all descendants of i, including i
    desc = c(i, descendants(x, i, internal=TRUE))
    A[fast.grid(rep(list(desc), 2))] = TRUE
  }
  A
}


#' @rdname ped_utils
#' @export
subnucs = function(x) {
  n = pedsize(x)
  if (is.singleton(x))
    return(list())
  FIDX = x$FIDX; MIDX = x$MIDX

  # Indices of unique parent couples
  p_pairs = paste(FIDX, MIDX)
  p_pairs_idx = which(!duplicated(p_pairs) & p_pairs != "0 0")

  # List all nucs: Format = c(father, mother, children)
  lapply(rev(p_pairs_idx), function(j) {
    nuc = list(father=FIDX[j], mother=MIDX[j], children=which(FIDX == FIDX[j] & MIDX == MIDX[j]))
    class(nuc) = "nucleus"
    attr(nuc, 'labels') = labels(x)
    nuc
  })
}


#' @rdname ped_utils
#' @export
peelingOrder = function(x) {
  # output: list of nuclear subfamilies. Format for each nuc:
  # list(father,mother,children,link), where link=0 for the last nuc.
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
      links = nucmembers[nucmembers %in% unlist(nucs[-i], use.names=F)]
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
    linktext = "0  (end of chain)"
  else {
    who = if (link == fa) "father"
    else if (link == mo) "mother"
    else if (link %in% ch) "child"
    else stop2("Erroneous link individual (not a member of nucleus): ", link)
    linktext = sprintf("%s (%s)", labs[link], who)
  }
  cat(sprintf("Pedigree nucleus.\nFather: %s\nMother: %s\nChildren: %s\nLink individual: %s\n",
              labs[fa], labs[mo], toString(labs[ch]), linktext))
}


has_numlabs = function(x) {
  # Returns TRUE if the labels of x are coercible to integers
  labs = labels(x)
  numlabs = suppressWarnings(as.character(as.integer(labs)))
  isTRUE(all(labs == numlabs))
}


.generations = function(x) {
  FOU = founders(x, internal=T)
  max(lengths(unlist(.descentPaths(x, FOU, internal=T), recursive=F)))
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
  nextNNnum = max(NNnum, na.rm=T) + 1
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


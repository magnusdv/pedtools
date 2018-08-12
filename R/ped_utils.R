
#' @export
has_unbroken_loops = function(x) {
  isTRUE(x$UNBROKEN_LOOPS)
}

#' @export
has_inbred_founders = function(x) {
  finb = x$FOUNDER_INBREEDING
  !is.null(finb) && any(finb > 0)
}

#' @export
labels.ped = function(object, ...) {
  object$LABELS
}

#' Pairwise common ancestors
#'
#' Computes a matrix A whose entry `A[i,j]` is TRUE if pedigree members i and j have a common ancestor, and FALSE otherwise.
#'
#' @param x a [ped()] object.
#'
#' @examples
#'
#' x = fullSibMating(3)
#' A = hasCA(x)
#' stopifnot(A[1,1], !A[1,2], all(A[3:8, 3:8]))
#'
#' @export
hasCA = function(x) {
  N = pedsize(x)
  A = matrix(FALSE, ncol=N, nrow=N)
  FOU = founders(x, internal=TRUE)
  for(i in FOU) {
    # vector of all descendants of i, including i
    desc = c(i, descendants(x, i, internal=TRUE))
    A[fast.grid(rep(list(desc), 2))] = TRUE
  }
  A
}


# Checks whether the labels of a ped oject are coercible to integers
has_numlabs = function(x) {
  labs = labels(x)
  numlabs = suppressWarnings(as.character(as.integer(labs)))
  isTRUE(all(labs == numlabs))
}


#'
#'
#' Return the sex of specified pedigree members.
#' @param x A `ped` object.
#' @param ids A character vector (or coercible to one) containing ID labels.
#'
#' @return An integer vector of the same length as `ids`, with entries 0
#'   (unknown), 1 (male) or 2 (female).
#'
#' @examples
#' x = nuclearPed(1)
#' stopifnot(all(getSex(x) == c(1,2,1)))
#'
#' @export
getSex = function(x, ids = labels(x)) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  x$SEX[internalID(x, ids)]
}

#' Pedigree size
#'
#' Return the number of pedigree members.
#'
#' @param x a [`ped`] object.
#'
#' @return A positive integer.
#' @export
#'
#' @examples
#' x = nuclearPed(1)
#' stopifnot(pedsize(x)==3)
#'
pedsize = function(x) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  length(x$ID)
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


#' Nuclear subsets of a pedigree
#'
#' @param x A `ped` object.
#' @param ... Not used.
#'
#' @return A list of `nucleus` objects. Each nucleus is a list with elements `father`, `mother` and `children`.
#'
#' @examples
#' x = cousinsPed(0, child=TRUE)
#' subnucs(x)
#'
#' @export
subnucs = function(x) {
  n = pedsize(x)
  if (n == 1)
    return(list())
  FID = x$FID; MID = x$MID

  # Indices of unique parent couples
  p_pairs = paste(FID, MID)
  p_pairs_idx = which(!duplicated(p_pairs) & p_pairs != "0 0")

  # List all nucs: Format = c(father, mother, children)
  lapply(rev(p_pairs_idx), function(j) {
    nuc = list(father=FID[j], mother=MID[j], children=which(FID == FID[j] & MID == MID[j]))
    class(nuc) = "nucleus"
    attr(nuc, 'labels') = labels(x)
    nuc
  })
}

#' @export
peelingOrder = function(x) {
  # output: list of nuclear subfamilies. Format for each nuc:
  # list(father,mother,children,link), where link=0 for the last nuc.
  # TODO: move to pedprobr package
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

#' @rdname subnucs
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

#' Inbreeding coefficients of founders
#'
#' Functions to get or set inbreeding coefficients for the pedigree founders.
#'
#' @param x A `ped` object
#' @param ids Any subset of `founders(x)`. If `ids` is missing in
#'   `founder_inbreeding()`, it is set to `founders(x)`.
#' @param named A logical: If TRUE, the output vector is named with the ID
#'   labels
#'
#' @return For `founder_inbreeding`, a numeric vector of the same length as
#' `ids`, containing the founder inbreeding coefficients.
#'
#' For `founder_inbreeding<-` the updated `ped` object is returned.
#'
#' @examples
#' x = nuclearPed(1)
#' founder_inbreeding(x, ids = '1') = 1
#' founder_inbreeding(x, named = TRUE)
#'
#' # Setting all founders at once (replacement value is recycled)
#' founder_inbreeding(x, ids = founders(x)) = 0.5
#' founder_inbreeding(x, named = TRUE)
#'
#' # Alternative syntax, using a named vector
#' founder_inbreeding(x) = c('1'=0.1, '2'=0.2)
#' founder_inbreeding(x, named = TRUE)
#'
#' @export
founder_inbreeding = function(x, ids, named = FALSE) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")

  fou = founders(x)

  if (missing(ids))
    ids = fou
  else if(any(!ids %in% fou)) {
    internalID(x, ids) # quick hack to catch unknown labels
    stop2("Pedigree member is not a founder: ", setdiff(ids, fou))
  }

  finb = x$FOUNDER_INBREEDING

  if(is.null(finb))
    finb = rep(0, length(fou))

  res = finb[match(ids, fou)]
  if(named)
    names(res) = ids
  res
}

#' @param value A numeric of the same length as `ids`, entries in the interval
#'   `[0, 1]`. If the vector is named, then the names are interpreted as ID
#'   labels of the founders whose inbreeding coefficients should be set. In this
#'   case, the `ids` argument should not be used. (See examples.)
#'
#'
#' @rdname founder_inbreeding
#' @export
`founder_inbreeding<-` = function(x, ids, value) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(!is.numeric(value))
    stop2("Inbreeding coefficients must be numeric: ", value)

  illegal = value < 0 | value > 1
  if(any(illegal))
    stop2("Inbreeding coefficients must be in the interval [0, 1]: ", value[illegal])

  if(missing(ids) && is.null(names(value)))
    stop2("When argument `ids` is missing, the replacement vector must be named")
  if(!missing(ids) && !is.null(names(value)))
    stop2("When the replacement vector is named, the `ids` argument must be missing")

  if(missing(ids)) ids = names(value)

  if(length(value) == 1)
    value = rep_len(value, length(ids))
  else if(length(ids) != length(value))
    stop2("Replacement vector must have length 1 or length(ids)")

  if(anyDuplicated(ids) > 0)
    stop2("Duplicated ID label: ", anyDuplicated(ids))

  fou = founders(x)
  if(any(!ids %in% fou)) {
    internalID(x, ids) # quick hack to catch unknown labels
    stop2("Pedigree member is not a founder: ", setdiff(ids, fou))
  }

  finb = x$FOUNDER_INBREEDING
  if(is.null(finb))
    finb = rep(0, length(fou))

  finb[match(ids, fou)] = value
  x$FOUNDER_INBREEDING = finb
  x
}

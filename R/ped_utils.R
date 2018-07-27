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
  numlabs = suppressWarnings(as.character(as.integer(x$LABELS)))
  isTRUE(all(x$LABELS == numlabs))
}


#'
#'
#' Return the sex of specified pedigree members.
#' @param x A `ped` object.
#' @param labels A character vector (or coercible to one) containing ID labels.
#'
#' @return An integer vector of the same length as `labels`, with entries 0
#'   (unknown), 1 (male) or 2 (female).
#'
#' @examples
#' x = nuclearPed(1)
#' stopifnot(getSex(x, 1) == 1)
#'
#' @export
getSex = function(x, labels) {
  assert_that(is.ped(x))
  x$SEX[internalID(x, labels)]
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
        if (!internal) res = lapply(res, function(int_ids) x$LABELS[int_ids])
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
    attr(nuc, 'labels') = x$LABELS
    nuc
  })
}

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

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
  if(isFALSE(all(x$LABELS == 1:x$NIND)))
    stop("This is currently only implemented for pedigrees with ordering 1,2,...")
  A = matrix(F, ncol=x$NIND, nrow=x$NIND)
  for(i in x$FOUNDERS) {
    # vector of all descendants of i, including i
    desc = c(i, descendants(x,i,internal=TRUE))
    A[fast.grid(rep(list(desc), 2))] = T
  }
  A
}


# Checks whether the labels of a ped oject are coercible to integers
has_numlabs = function(x) {
  assert_that(is.ped(x))
  numlabs = suppressWarnings(as.numeric(x$LABELS))
  isTRUE(all(numlabs == as.integer(numlabs)))
}


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
#' stopifnot(pedSize(x)==3)
#'
pedSize = function(x) {
  assert_that(is.ped(x))
  length(x$ID)
}

catLabels = function(x, int_ids) paste(x$LABELS[int_ids], collapse = ", ")

.generations = function(x) {
    max(lengths(unlist(.descentPaths(x, x$FOUNDERS, internal = TRUE), recursive = F)))
}

.descentPaths = function(x, ids, internal = FALSE) {
    if (!internal) ids = internalID(x, ids)

    offs = lapply(1:pedSize(x), children, x = x, internal = TRUE)
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

#' @export
peeling_order = function(x) {
  # output: list of nuclear subfamilies. Format for each nuc:
  # list(father,mother,children,link), where link=0 for the last nuc.
  n = pedSize(x)
  if (n == 1)
    return(list())
  FID = x$FID; MID = x$MID

  # Indices of unique parent couples
  p_pairs = paste(FID, MID)
  p_pairs_idx = which(!duplicated(p_pairs) & p_pairs != "0 0")

  # List all nucs: Format = c(father, mother, children)
  list1 = lapply(rev(p_pairs_idx), function(j) {
    nuc = list(father=FID[j], mother=MID[j], children=which(FID == FID[j] & MID == MID[j]))
    class(nuc) = "nucleus"
    attr(nuc, 'labels') = x$LABELS
    nuc
  })

  peeling = vector("list", Nnucs <- length(list1))
  i = k = 1

  while (length(list1)) {
    # Start searching for a nuc with only one link to others
    nuc = list1[[i]]

    # Identify links to other remaining nucs
    if(length(list1) > 1) {
      nucmembers = c(nuc$father, nuc$mother, nuc$children)
      links = nucmembers[nucmembers %in% unlist(list1[-i], use.names=F)]
    }
    else
      links = 0 # if nuc is the last

    # If only one link: move nuc to peeling list, and proceed
    if (length(links) == 1) {
      nuc$link = links
      peeling[[k]] = nuc
      list1[i] = NULL
      i = 1
      k = k+1
    }
    else {
      if (i == length(list1)) { # LOOP! Include remaining nucs without 'link', and break.
        peeling[k:Nnucs] = list1
        break
      }
      # Otherwise try next remaining nuc
      i = i+1
    }
  }

  peeling
}

#' @export
print.nucleus = function(nuc) {
  labs = attr(nuc, 'labels')
  fa = nuc$father
  mo = nuc$mother
  ch = nuc$children

  link = nuc$link
  if (is.null(link))
    linktext = "PART OF LOOP; NO LINK ASSIGNED"
  else if (link==0)
    linktext = "NONE - end of chain"
  else {
    who = if (link == fa) "father"
          else if (link == mo) "mother"
          else if (link %in% ch) "child"
          else stop("link individual ", link, " is not a member of the nucleus")
    linktext = sprintf("%d (%s)", link, who)
  }
  cat(sprintf("Pedigree nucleus.\nFather: %s\nMother: %s\nChildren: %s\nLink individual: %s\n",
              labs[fa], labs[mo], paste(labs[ch], collapse=", "), linktext))
}

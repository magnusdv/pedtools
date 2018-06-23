#' Pairwise common ancestors
#'
#' Computes a matrix A whose entry A[i,j] is TRUE if pedigree members i and j have a common ancestor, and FALSE otherwise.
#'
#' @param x a \code{\link{ped}} object.
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
  numlabs = suppressWarnings(as.numeric(x$LABELS))
  isTRUE(all(numlabs == as.integer(numlabs)))
}


getSex = function(x, labels)
  x$SEX[internalID(x, labels)]

pedSize = function(x) length(x$ID)

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


peeling_order = function(x) {
  # output: list of nuclear subfamilies. Format for each nuc:
  # c(link,father,mother,offsp1,..), where pivot=0 for the last nuc.
  n = pedSize(x)
  if (n == 1)
    return(list())
  FID = x$FID; MID = x$MID

  # Indices of unique parent couples
  p_pairs = paste(FID, MID)
  p_pairs_idx = which(!duplicated(p_pairs) & p_pairs != "0 0")

  # List all nucs: Format = c(father, mother, children)
  list1 = lapply(rev(p_pairs_idx), function(j) {
    nuc = c(father=FID[j], mother=MID[j], child=which(FID == FID[j] & MID == MID[j]))
    class(nuc) = "nucleus"
    nuc})

  peeling = vector("list", Nnucs <- length(list1))
  i = k = 1

  while (length(list1)) {
    # Start searching for a nuc with only one link to others
    nuc = list1[[i]]

    # Identify links to other remaining nucs
    if(length(list1) > 1)
      links = nuc[nuc %in% unlist(list1[-i])]
    else
      links = 0 # if nuc is the last

    # If only one link: move nuc to peeling list, and proceed
    if (length(links) == 1) {
      attr(nuc, 'link') = links
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

print.nucleus = function(nuc) {
  children = paste(nuc[-(1:2)], collapse=", ")
  link = attr(nuc, 'link')
  linktext = if (is.null(link)) "PART OF LOOP; NO LINK ASSIGNED"
  else if (link==0) "End of chain"
  else sprintf("%d (%s)", link, switch(link, "father", "mother", "child"))
  cat(sprintf("Father: %d\nMother: %d\nChildren: %s\nLink: %s\n",
              nuc[['father']], nuc[['mother']], children, linktext))
}

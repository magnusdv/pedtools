#' Random pedigree
#'
#' Generate a random connected pedigree by applying random mating starting from
#' a finite population.
#'
#' Starting from an initial set of `f` singletons, a sequence of `n-f` random
#' matings is performed. The sampling of parents in each mating is set up to
#' ensure that the final result is connected.
#'
#' @param n A positive integer: the total number of individuals. Must be at
#'   least 3.
#' @param f A positive integer: the number of founders. Must be at least 2
#'   unless selfing is allowed.
#' @param maxDirectGap An integer; the maximum distance between direct
#'   descendants allowed to mate. For example, the default value of 1 allows
#'   parent-child mating, but not grandparent-grandchild. Set to `Inf` or `NULL`
#'   to remove all restrictions.
#' @param selfing A logical indicating if selfing is allowed. Default: FALSE.
#' @param seed An integer seed for the random number generator (optional).
#' @param g,founders Deprecated arguments.
#'
#' @return A connected pedigree returned as a `ped` object.
#'
#' @examples
#' plot(randomPed(7, f = 2, seed = 12))
#'
#' # Disallow mating between direct discendants
#' plot(randomPed(7, f = 2, seed = 12, maxDirectGap = 0))
#'
#' # No restrictions on mating between direct discendants
#' plot(randomPed(7, f = 2, seed = 12, maxDirectGap = Inf))
#'
#' # Allow selfing
#' y = randomPed(5, f = 2, seed = 2, selfing = TRUE)
#' hasSelfing(y)
#' y
#' plot(y, arrows = TRUE)
#'
#' @export
randomPed = function(n, f = 2, maxDirectGap = 1, selfing = FALSE,
                     seed = NULL, g = NULL, founders = NULL) {
  # TODO: Remove
  if(!is.null(g)) {
    message("Switching to the old version, using deprecated argument `g`.")
    return(.randomPed(g, founders %||% f, selfing = selfing, seed = seed))
  }

  if(!is.null(seed))
    set.seed(seed)

  if(!isCount(n))
    stop2("Argument `n` must be a natural number: ", n)
  if(n < 3)
    stop2("The total number of individuals must be at least 3: ", n)
  if(!isCount(f))
    stop2("Argument `f` must be a natural number")
  if(f < 2 && !selfing)
    stop2("When selfing is disallowed, the number of founders must be at least 2: ", f)
  if(f > (n+1)/2)
    stop2(sprintf("Too many founders; the pedigree cannot be connected.\nNote: A pedigree with %d members can have at most %d founders.",
                  n, floor((n+1)/2)))

  id = seq_len(n)
  fid = mid = numeric(n)

  # Sex of founders
  if(selfing && f==1) {
    sexFou = sample.int(2, 1)
  }
  else {
    sexFou = rep(2L, f)
    nMaleFou = sample.int(f-1, 1)
    sexFou[sample.int(f, size = nMaleFou)] = 1L
  }

  # Sex of nonfounders
  sexNonfou = sample.int(2, size = n - f, replace = TRUE)

  # Complete sex vector
  sex = c(sexFou, sexNonfou)

  # Component index.
  # Initially 1,2,...f (singletons) and unknown (0) for the rest
  compvec = c(seq_len(f), rep(0L, n-f))

  # Generation gap matrix
  gengap = matrix(NA_integer_, nrow = n, ncol = n)

  for(k in seq(f + 1, to = n)) {

    sq = seq_len(k-1)
    ncomp = max(compvec)
    compSizes = tabulate(compvec)

    # Minimum number of pairings required to make output connected
    requ = ncomp-1

    # Surplus pairings (i.e., leeway in addition to the required ones)
    surplus = (n-k+1) - requ

    # Weights for choosing the first parent
    # The fewer surplus pairings, the more weight on small components.
    # No surplus: Must take smallest (typically singleton)
    wComp = (min(compSizes)/compSizes)^(1/surplus)
    w1 = wComp[compvec[sq]]
    w1 = w1/sum(w1)

    # Candidates to avoid for first parent (maxDirectGap may be prohibitive)
    if(!selfing && !is.null(maxDirectGap) && maxDirectGap < Inf) {
       avoid = sapply(sq, function(i) {
         gapsi = pmax(gengap[sq, i], gengap[i, sq], na.rm = TRUE)
         toofar = !is.na(gapsi) & gapsi > maxDirectGap
         samesex = sex[sq] == sex[i]
         all(toofar | samesex)
       })
    }
    else {
      avoid = rep(FALSE, k-1)
    }

    # Sample first parent!
    cand1 = sq[!avoid]
    prob1 = w1[!avoid]
    par1 = sample(cand1, size = 1, prob = prob1)
    comp1 = compvec[par1]

    # Candidates for the other parent: Any of opposite sex
    isCand = sex[sq] != sex[par1]

    # If selfing, include self as candidate
    if(selfing)
      isCand[par1] = TRUE

    # If no surplus, partner cannot be in same component
    if(surplus == 0)
      isCand[compvec[sq] == comp1] = FALSE

    # Apply generation gap limit if given
    if(!is.null(maxDirectGap) && maxDirectGap < Inf) {
      gaps = pmax(gengap[sq, par1], gengap[par1, sq], na.rm = TRUE)
      toofar = !is.na(gaps) & gaps > maxDirectGap
      isCand[toofar] = FALSE
    }

    # If no candidates, plot (for debugging)
    if(!any(isCand)) {
      p = ped(id[sq], fid[sq], mid[sq], sex[sq])
      plot(p, hatched = par1)
    }

    # The fewer surplus, the more probable to go to another comp
    par2 = safe_sample(sq[isCand], 1)
    comp2 = compvec[par2]

    # Swap if par1 female & par2 male
    if(sex[par1] > sex[par2]) {
      fid[k] = par2
      mid[k] = par1
    } else {
      fid[k] = par1
      mid[k] = par2
    }

    # If selfing: set sex to 0
    if(par1 == par2)
      sex[k] = 0L

    # Modify components
    if(comp1 == comp2)
      compvec[k] = comp1
    else {
      compvec[k] = cmp = min(comp1, comp2)
      compvec[compvec == max(comp1, comp2)] = cmp

      # Reduce larger comp numbers by 1
      b = compvec > max(comp1, comp2)
      compvec[b] = compvec[b] - 1
    }

    # Update generation gap matrix
    gengap[sq, k] = 1L + pmax(gengap[sq, par1], gengap[sq, par2], na.rm = TRUE)
    if(is.na(gengap[par1, k]))
      gengap[par1, k] = 1L
    if(is.na(gengap[par2, k]))
      gengap[par2, k] = 1L
  }

  newPed(as.character(id), as.integer(fid), as.integer(mid), sex, FAMID = "")
}

#' @importFrom stats rpois
.randomPed = function(g, founders = rpois(1,3)+1, selfing = FALSE, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  if(!selfing) founders = max(2, founders)
  id = seq_len(founders + g)
  fid = mid = numeric(founders + g)


  foundersM = ceiling(founders/2)
  foundersF = founders - foundersM
  sex = c(rep(1:2, c(foundersM, foundersF)), sample.int(2, size = g, replace = TRUE))
  males = (sex == 1)
  females = (sex == 2)

  for(k in seq(founders+1, length = g)) {

    if(selfing) {
      # First parent: Any member
      par1 = sample.int(k-1, 1)

      # Second parent: Either par1 or any of opposite sex. Note: k > 1.
      par2_candidates = c(par1, which(sex[1:(k-1)] != sex[par1]))
      par2 = safe_sample(par2_candidates, 1)

      # Swap if par1 female & par2 male
      if(sex[par1] > sex[par2]) {
        fid[k] = par2
        mid[k] = par1
      } else {
        fid[k] = par1
        mid[k] = par2
      }
      if(par1 == par2) sex[k] = 0
    } else {
      potential_fathers = which(males[1:(k-1)])
      potential_mothers = which(females[1:(k-1)])
      fid[k] = safe_sample(potential_fathers, 1)
      mid[k] = safe_sample(potential_mothers, 1)
    }
  }
  ped(id, fid, mid, sex, validate = FALSE, reorder = FALSE)
}

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
#' @param selfing A logical indicating if selfing is allowed. Default: FALSE.
#' @param seed An integer seed for the random number generator (optional).
#' @param g,founders Deprecated arguments.
#'
#' @return A connected pedigree returned as a `ped` object.
#'
#' @examples
#' randomPed(8, f = 3, seed = 11)
#' randomPed(8, f = 3, seed = 11, selfing = TRUE)
#'
#' @export
randomPed = function(n, f = 2, selfing = FALSE, seed = NULL, g = NULL, founders = NULL) {
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

  # List of components. Initially f singletons
  compvec = c(seq_len(f), rep(0L, n-f))

  for(k in seq(f + 1, to = n)) {

    sq = seq_len(k-1)
    ncomp = max(compvec)
    compSizes = tabulate(compvec)

    # Surplus pairings
    surplus = (n-k+1) - (ncomp-1)

    # Weights for choosing the first parent
    # The fewer surplus pairings, the more weight on small components.
    # No surplus: Must take smallest (typically singleton)
    wComp = (min(compSizes)/compSizes)^(1/surplus)
    w1 = wComp[compvec[sq]]
    w1 = w1/sum(w1)

    # Sample first parent!
    par1 = sample.int(k-1, 1, prob = w1)
    comp1 = compvec[par1]

    # Second parent: Any of opposite sex.
    isCand = sex[sq] != sex[par1]
    if(selfing)
      isCand[par1] = TRUE

    # If no surplus, partner cannot be in same component
    if(surplus == 0)
      isCand[compvec[sq] == comp1] = FALSE

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

#' Random pedigree
#'
#' Generate a random pedigree by applying random mating starting from a finite
#' population. The resulting pedigree will have `f + g` members, where `f` is
#' the number of founders and `g` is the number of matings.
#'
#' The sampling scheme for choosing parents in each mating depends on the
#' `selfing` parameter. If `selfing = FALSE`, a father is randomly sampled from
#' the existing males, and a mother from the existing females. If `selfing =
#' TRUE` then one parent P1 is sampled first (among all members), and then a
#' second parent from the set consisting of P1 and all members of the opposite
#' sex. The gender of the child is randomly chosen with equal probabilities.
#'
#' @param g A positive integer: The number of matings.
#' @param founders A positive integer: The size of the initial population.
#' @param selfing A logical indicating if selfing is allowed.
#' @param seed A numerical seed for random number generation. (Optional.)
#' @return A `ped` object.
#'
#' @examples
#' randomPed(3, 3)
#' randomPed(3, 3, selfing = TRUE)
#'
#' @importFrom stats rpois
#' @export
randomPed = function(g, founders = rpois(1,3)+1, selfing = FALSE, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  if(!selfing) founders = max(2, founders)
  id = seq_len(founders + g)
  fid = mid = numeric(founders + g)

  foundersM = ceiling(founders/2)
  foundersF = founders - foundersM
  sex = c(rep(1:2, c(foundersM, foundersF)), sample.int(2, size = g, replace = T))
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

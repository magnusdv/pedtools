safe_sample <- function(x, ...) x[sample.int(length(x), ...)]

#' Random pedigree
#'
#' Generate a random pedigree by applying random mating starting from a finite
#' population. The resulting pedigree will have \code{f + g} members, where
#' \code{f} is the number of founders and \code{g} is the number of matings.
#'
#' The sampling scheme for chosing parents in each mating depends on the
#' \code{selfing} parameter. If \code{selfing=FALSE}, a father is randomly
#' sampled from the exisiting males, and a mother from the existing females. If
#' \code{selfing=TRUE} then one parent P1 is sampled first (among all members),
#' and then a second parent from the set consisting of P1 and all members of the
#' opposite sex. The gender of the child is randomly chosen with equal
#' probabilities.
#'
#' @param g A positive integer: The number of matings.
#' @param founders A positive integer: The size of the initial population.
#' @param selfing A logical indicating if selfing is allowed.
#' @param seed A numerical seed for random number generation. (Optional.)
#' @return A \code{ped} object.
#'
#' @examples
#' randomPed(3, 3)
#' randomPed(3, 3, selfing=TRUE)
#'
#' @importFrom stats rpois
#' @export
randomPed = function(g, founders=rpois(1,3)+1, selfing=FALSE, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  if(!selfing) founders = max(2, founders)
  ID = seq_len(founders + g)
  FID = MID = numeric(founders + g)

  foundersM = ceiling(founders/2)
  foundersF = founders - foundersM
  SEX = c(rep(1:2, c(foundersM, foundersF)), sample.int(2, size=g, replace = T))
  males = (SEX == 1)
  females = (SEX == 2)

  for(k in seq(founders+1, length=g)) {

    if(selfing) {
      # First parent: Any member
      par1 = sample.int(k-1, 1)

      # Second parent: Either par1 or any of opposite sex. Note: k > 1.
      par2_candidates = c(par1, which(SEX[1:(k-1)] != SEX[par1]))
      par2 = safe_sample(par2_candidates, 1)

      # Swap if par1 female & par2 male
      if(SEX[par1] > SEX[par2]) {
        FID[k] = par2
        MID[k] = par1
      } else {
      FID[k] = par1
      MID[k] = par2
      }
      if(par1==par2) SEX[k] = 0
    } else {
      potential_fathers = which(males[1:(k-1)])
      potential_mothers = which(females[1:(k-1)])
      FID[k] = safe_sample(potential_fathers, 1)
      MID[k] = safe_sample(potential_mothers, 1)
    }
  }
  ped(ID, FID, MID, SEX, check=FALSE, reorder=FALSE)
}

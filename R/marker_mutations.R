#' Mutation matrices
#'
#' Construct mutation matrices
#'
#' At the moment only `model = equal` is implemented
#' @param alleles A character vector (or coercible to character) with allele labels
#' @param model A string: either "equal", "proportional" or "random"
#' @param rate A number between 0 and 1
#'
#' @return A square matrix whose `dimnames` are given by the `alleles`.
#'
#' @examples
#' mutationMatrix(alleles = 1:3, model = "equal", rate = 0.05)
#'
#' @export
mutationMatrix = function(alleles, model = c("equal", "proportional", "random"), rate) {
  nall = length(alleles)
  if(nall < 2)
    return(NULL)
  mutmat = matrix(ncol = nall, nrow = nall, dimnames = list(alleles, alleles))

  model = match.arg(model)
  if(model == "equal") {
    mutmat[] = rate/(nall - 1)
    diag(mutmat) = 1 - rate
  }
  else
    stop("Only `model=equal` is implemented")
  mutmat
}

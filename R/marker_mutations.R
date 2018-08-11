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


.checkMutationMatrix = function(mutmat, alleles) {
  # Check that mutation matrix is compatible with allele number / names.  Sort matrix
  # according to the given allele order (this is important since dimnames are NOT used in
  # calculations).
  alleles = as.character(alleles)
  N = length(alleles)

  if (any((dm <- dim(mutmat)) != N)) {
    print(mutmat)
    stop2(sprintf("Dimension of mutation matrix (%d x %d)
                  inconsistent with number of alleles (%d)", dm[1], dm[2], N))
  }

  if (any(round(.rowSums(mutmat, N, N), 3) != 1)) {
    print(mutmat)
    stop2("Row sums of mutation matrix are not 1")
  }

  row_nms = rownames(mutmat)
  col_nms = colnames(mutmat)
  if (!setequal(row_nms, alleles) || !setequal(col_nms, alleles)) {
    print(mutmat)
    stop2("Dimnames of mutation do not match allele names")
  }

  m = mutmat[alleles, alleles]

  # Lumbability: always lumpable (any alleles) if rows are identical (except diagonal)
  attr(m, "lumpability") = NA
  if (N > 2) {
    if (all(vapply(1:N, function(i) diff(range(m[-i, i])) == 0, logical(1))))
      attr(m, "lumpability") = "always"  # If each column has 0 range
  }
  m
  }

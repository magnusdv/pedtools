#' Mutation matrices
#'
#' Construct mutation matrices.
#'
#' Descriptions of the models:
#'
#' * "equal" :  All mutations equally likely; probability `1-rate` of no mutation.
#' * "proportional" : Mutation probabilities are proportional to the target allele frequencies
#' * "random" : This produces a matrix of random numbers, where each row is normalised so that it sums to 1
#'
#' @param alleles A character vector (or coercible to character) with allele labels
#' @param model A string: either "equal", "proportional" or "random"
#' @param afreq A numeric vector of allele frequencies. Required in model "proportional"
#' @param rate A number between 0 and 1. Required in models "equal" and "proportional"
#' @param seed A single number. Optional parameter in the "random" model, passed on to `set.seed()`
#'
#' @return A square matrix with entries in `[0, 1]`, with dimnames taken from `alleles`.
#'
#' @examples
#' mutationMatrix(alleles = 1:3, model = "equal", rate = 0.05)
#'
#' @importFrom stats runif
#' @export
mutationMatrix = function(alleles, model = c("equal", "proportional", "random"),
                          afreq = NULL, rate = NULL, seed = NULL) {
  nall = length(alleles)
  if(nall < 2)
    return(NULL)
  mutmat = matrix(ncol = nall, nrow = nall, dimnames = list(alleles, alleles))

  switch(match.arg(model),
    equal = {
      mutmat[] = rate/(nall - 1)
      diag(mutmat) = 1 - rate
    },
    proportional = {
      alpha = rate / sum(afreq * (1 - afreq))
      mutmat[] = (1 - alpha) * diag(nall) + alpha * rep(afreq, each = nall)
      if(max(mutmat) > 1 || min(mutmat) < 0)
        stop("Impossible mutation matrix. Reduce rate")
    },
    random = {
      if(!is.null(seed))
        set.seed(seed)
      mutmat[] = runif(nall^2, min = 0, max = 1)
      mutmat = mutmat / rowSums(mutmat)
    })
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


#' Create simple pedigrees
#'
#' These are utility functions for creating some common pedigree structures as
#' `ped` objects.
#'
#' `halfSibPed(nch1, nch2)` produces a pedigree containing two sibships (of
#' sizes `nch1` and `nch2`) with the same father, but different mothers. If
#' maternal halfsibs are wanted instead, use [swapSex()] afterwards. (See
#' examples below.)
#'
#' The call `cousinsPed(degree=n, removal=k)` creates a pedigree with two n'th
#' cousins, k times removed. By default, removals are added on the right side.
#' To override this, the parameter `degree2` can be used to indicate explicitly
#' the number of generations on the right side of the pedigree. When `degree2`
#' is given `removal` is ignored. (Similarly for `halfCousinsPed`.)
#'
#' @param nch The number of children. If NULL, it is taken to be the
#'   `length(children)`
#' @param sex A vector with integer gender codes (0=unknown, 1=male, 2=female).
#'   In `nuclearPed()`, it contains the genders of the children and is recycled
#'   (if neccessary) to length `nch`. In `linearPed()` it also contains the
#'   genders of the children (1 in each generation) and should have length at
#'   most `n` (recycled if shorter than this).
#' @param father The label of the father.
#' @param mother The label of the father.
#' @param children A character of length `nch`, with labels of the children.
#' @param nch1,nch2 The number of children in each sibship.
#' @param sex1,sex2 Vectors of gender codes for the children in each sibship.
#'   Recycled (if neccessary) to lengths `nch1` and `nch2` respectively.
#' @param n The number of generations, not including the initial founders.
#' @param degree,degree2 Non-negative integers, indicating the degree of
#'   cousin-like relationships: 0=siblings, 1=first cousins; 2=second cousins,
#'   a.s.o. See Details and Examples.
#' @param removal Non-negative integers, indicating the removals of cousin-like
#'   relationships. See Details and Examples.
#' @param child A logical: Should an inbred child be added to the two cousins?
#'
#' @return A `ped` object.
#'
#' @seealso [ped()], [singleton()], [ped_complex], [ped_subgroups]
#'
#' @examples
#'
#' # A nuclear family with 2 boys and 3 girls
#' nuclearPed(5, sex=c(1,1,2,2,2))
#'
#' # A straight line of females
#' linearPed(3, sex = 2)
#'
#' # Paternal half brothers
#' x = halfSibPed()
#'
#' # Change into maternal half brothers
#' x = swapSex(x, 1)
#'
#' # Larger half sibships: boy and girl on one side; 3 girls on the other
#' halfSibPed(nch1 = 2, sex = 1:2, nch2 = 3, sex2 = 2)
#'
#' # Grand aunt:
#' cousinsPed(degree=0, removal=2)
#'
#' # Second cousins once removed.
#' cousinsPed(degree=2, removal=1)
#'
#' # Same, but with the 'removal' on the left side.
#' cousinsPed(degree=3, degree2=2)
#'
#' # A child of first cousin parents.
#' cousinsPed(degree=1, child=TRUE)
#'
#' @name ped_basic
NULL

#' @rdname ped_basic
#' @export
nuclearPed = function(nch, sex = 1, father = '1', mother = '2',
                      children = as.character(seq.int(3, length.out=nch))) {
  if(missing(nch))
    nch = length(children)
  if(!is_count(nch))
    stop2("`nch` must be a positive integer: ", nch)
  if(length(children) != nch)
    stop2("`children` must have length `nch`")
  if(length(father) != 1)
    stop2("`father` must have length 1")
  if(length(mother) != 1)
    stop2("`mother` must have length 1")
  if(father == "1" && "1" %in% c(mother, children))
    stop2("By default the father is named '1'. ",
          "If you want to use this for someone else, please specify a different label for the father.")
  if(mother == "2" && "2" %in% c(father, children))
    stop2("By default the mother is named '2'. ",
          "If you want to use this for someone else, please specify a different label for the mother.")

  sex = validate_sex(sex, nInd = nch)

  x = ped(id = 1:(2 + nch),
      fid = c(0, 0, rep(1, nch)),
      mid = c(0, 0, rep(2, nch)),
      sex = c(1, 2, sex))

  relabel(x, c(father, mother, children))
}

#' @rdname ped_basic
#' @export
halfSibPed = function(nch1 = 1, nch2 = 1, sex1 = 1, sex2 = 1) {
  if(!is_count(nch1))
    stop2("`nch1` must be a positive integer: ", nch1)
  if(!is_count(nch2))
    stop2("`nch2` must be a positive integer: ", nch2)
  sex1 = validate_sex(sex1, nInd = nch1)
  sex2 = validate_sex(sex2, nInd = nch2)

  x = nuclearPed(nch1, sex = sex1)
  x = relabel(x, c(1,2, 4:(1 + pedsize(x))))
  x = addChildren(x, father = 1, mother = 3, nch = nch2, sex = sex2,
                  verbose = F)
  x
}

#' @rdname ped_basic
#' @export
linearPed = function(n, sex = 1) {
  if(!is_count(n))
    stop2("`n` must be a positive integer: ", n)

  sex = validate_sex(sex, nInd = n)

  nInd = 1 + 2*n
  child_idx = seq(3, nInd, by=2)

  # Create ped
  id = 1:nInd
  fid = mid = rep_len(0, nInd)
  fid[child_idx] = 2*(1:n) - 1
  mid[child_idx] = 2*(1:n)
  sex0 = rep_len(1:2, nInd)

  x = ped(id, fid, mid, sex0, validate = F,verbose = F, reorder = F)

  # swap genders if needed
  if(any(sex == 2)) {
    swaps = child_idx[sex == 2]
    x = swapSex(x, swaps, verbose = FALSE)
  }

  x
}

#' @rdname ped_basic
#' @export
cousinsPed = function(degree, removal = 0, degree2 = NULL, child = FALSE) {
  # Creates a pedigree linking two n'th cousins k times removed, where
  # n=degree, k=removal. By default, removals are added on the right side,
  # i.e. degree2 = degree + removal.  If degree2 is non-NULL, the removal
  # parameter is ignored.
  if(!is_count(degree, minimum=0)) stop2("`degree` must be a nonnegative integer: ", degree)
  if(!is_count(removal, minimum=0)) stop2("`removal` must be a nonnegative integer: ", removal)
  if(!is.null(degree2) && !is_count(degree2, minimum=0))
    stop2("`degree` must be a positive integer (or NULL): ", degree2)

  if (is.null(degree2))
      degree2 = degree + removal

  # Chain on the left side
  x = nuclearPed(1)
  for (i in seq_len(degree)) x = addSon(x, pedsize(x), verbose = F)

  # Chain on the right side
  x = addChildren(x, father = 1, mother = 2, nch = 1, verbose = F)
  for (i in seq_len(degree2)) x = addSon(x, pedsize(x), verbose = F)
  x = swapSex(x, pedsize(x), verbose = F)

  if (child) {
      cous = leaves(x)
      x = addChildren(x, father = cous[1], mother = cous[2], nch = 1, verbose = F)
  }
  x
}

#' @rdname ped_basic
#' @export
halfCousinsPed = function(degree, removal = 0, degree2 = NULL, child = FALSE) {
  # Creates a pedigree linking two n'th half cousins k times removed, where
  # n=degree, k=removal.  By default, removals are added on the right side,
  # i.e. degree2 = degree + removal.  If degree2 is non-NULL, the removal
  # parameter is ignored.
  if(!is_count(degree, minimum=0)) stop2("`degree` must be a nonnegative integer: ", degree)
  if(!is_count(removal, minimum=0)) stop2("`removal` must be a nonnegative integer: ", removal)
  if(!is.null(degree2) && !is_count(degree2, minimum=0))
    stop2("`degree` must be a positive integer (or NULL): ", degree2)

  if (is.null(degree2))
    degree2 = degree + removal

  # Chain on the left side
  x = nuclearPed(1)
  for (i in seq_len(degree)) x = addSon(x, pedsize(x), verbose = F)

  # Chain on the right side
  x = addSon(x, 1, verbose = F)
  for (i in seq_len(degree2)) x = addSon(x, pedsize(x), verbose = F)
  x = swapSex(x, pedsize(x), verbose = F)

  if (child) {
    cous = leaves(x)
    x = addChildren(x, father = cous[1], mother = cous[2], nch = 1, verbose = F)
  }
  x
}


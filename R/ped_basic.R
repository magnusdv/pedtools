#' Create simple pedigrees
#'
#' Utility functions for creating some common pedigree structures.
#'
#' `halfSibPed(nch1, nch2)` produces a pedigree containing two sibships (of
#' sizes `nch1` and `nch2`) with the same father, but different mothers. If
#' maternal halfsibs are wanted instead, use [swapSex()] afterwards. (See
#' examples below.)
#'
#' `cousinPed(degree = n, removal = k)` creates a pedigree with two `n`'th
#' cousins, `k` times removed. By default, removals are added on the right side,
#' but this can be changed by adding `side = left`. (Similarly for
#' `halfCousinPed`.)
#'
#' `ancestralPed(g)` returns the family tree of a single individual, including
#' all ancestors `g` generations back.
#'
#' `selfingPed(s)` returns a line of `s` consecutive selfings.

#' @param nch The number of children. If NULL, it is taken to be the
#'   `length(children)`
#' @param sex A vector with integer gender codes (0=unknown, 1=male, 2=female).
#'   In `nuclearPed()`, it contains the genders of the children and is recycled
#'   (if necessary) to length `nch`. In `linearPed()` it also contains the
#'   genders of the children (1 in each generation) and should have length at
#'   most `n` (recycled if shorter than this). In `selfingPed()` it should be a
#'   single number, indicating the gender of the last individual (the others
#'   must necessarily have gender code 0.)
#' @param father The label of the father.
#' @param mother The label of the mother.
#' @param children A character of length `nch`, with labels of the children.
#' @param nch1,nch2 The number of children in each sibship.
#' @param sex1,sex2 Vectors of gender codes for the children in each sibship.
#'   Recycled (if necessary) to lengths `nch1` and `nch2` respectively.
#' @param n The number of generations, not including the initial founders.
#' @param degree A non-negative integer: 0=siblings, 1=first cousins; 2=second
#'   cousins, a.s.o.
#' @param removal A non-negative integer. See Details and Examples.
#' @param side Either "right" or "left"; the side on which removals should be
#'   added.
#' @param child A logical: Should an inbred child be added to the two cousins?
#' @param g A nonnegative integer indicating the number of ancestral generations
#'   to include. The resulting pedigree has `2^(g+1)-1` members. The case `g =
#'   0` results in a singleton.
#' @param s A nonnegative integer indicating the number of consecutive selfings.
#'   The case `s = 0` results in a singleton.
#' @return A `ped` object.
#'
#' @seealso [ped()], [singleton()], [ped_complex], [ped_subgroups]
#'
#' @examples
#'
#' # A nuclear family with 2 boys and 3 girls
#' nuclearPed(5, sex = c(1, 1, 2, 2, 2))
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
#' cousinPed(degree = 0, removal = 2)
#'
#' # Second cousins once removed.
#' cousinPed(degree = 2, removal = 1)
#'
#' # Same, but with the 'removal' on the left side.
#' cousinPed(2, 1, side = "left")
#'
#' # A child of half first cousins.
#' halfCousinPed(degree = 1, child = TRUE)
#'
#' # The 'family tree' of a person
#' ancestralPed(g = 2)
#'
#' @name ped_basic
NULL

#' @rdname ped_basic
#' @export
nuclearPed = function(nch, sex = 1, father = '1', mother = '2',
                      children = as.character(seq.int(3, length.out = nch))) {
  if(missing(nch))
    nch = length(children)
  if(!isCount(nch))
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
  if(!isCount(nch1))
    stop2("`nch1` must be a positive integer: ", nch1)
  if(!isCount(nch2))
    stop2("`nch2` must be a positive integer: ", nch2)
  sex1 = validate_sex(sex1, nInd = nch1)
  sex2 = validate_sex(sex2, nInd = nch2)

  x = ped(id = seq_len(3 + nch1 + nch2),
          fid = c(0, 0, 0, rep.int(1, nch1 + nch2)),
          mid = c(0, 0, 0, rep.int(2, nch1), rep.int(3, nch2)),
          sex = c(1, 2, 2, sex1, sex2),
          verbose = FALSE)
  x
}

#' @rdname ped_basic
#' @export
linearPed = function(n, sex = 1) {
  if(!isCount(n, minimum = 0))
    stop2("`n` must be a nonnegative integer: ", n)

  sex = validate_sex(sex, nInd = n)

  nInd = 1 + 2*n
  child_idx = seq(3, nInd, by = 2)

  # Create ped
  id = 1:nInd
  fid = mid = rep_len(0, nInd)
  fid[child_idx] = 2*(1:n) - 1
  mid[child_idx] = 2*(1:n)
  sex0 = rep_len(1:2, nInd)

  x = ped(id, fid, mid, sex0, validate = FALSE,verbose = FALSE, reorder = FALSE)

  # swap genders if needed
  if(any(sex == 2)) {
    swaps = child_idx[sex == 2]
    x = swapSex(x, swaps, verbose = FALSE)
  }

  x
}

#' @rdname ped_basic
#' @export
cousinPed = function(degree, removal = 0, side = c("right", "left"), child = FALSE) {
  if(!isCount(degree, minimum = 0))
    stop2("`degree` must be a nonnegative integer: ", degree)
  if(!isCount(removal, minimum = 0))
    stop2("`removal` must be a nonnegative integer: ", removal)

  deg_right = deg_left = degree
  switch(match.arg(side),
         right = {deg_right <- deg_right + removal},
         left  = {deg_left <- deg_left + removal})

  # Chain on the left side
  x = linearPed(deg_left + 1)

  # Chain on the right side
  y = linearPed(deg_right + 1)
  y = relabel(y, old = 3:pedsize(y),
              new = pedsize(x) + 1:(pedsize(y) - 2))

  # Merge
  z = mergePed(x, y)

  if (child) {
    parents = leaves(z)
    z = swapSex(z, parents[2], verbose = FALSE)
    z = addChildren(z, father = parents[1], mother = parents[2], nch = 1, verbose = FALSE)
  }
  z
}


#' @rdname ped_basic
#' @export
halfCousinPed = function(degree, removal = 0, side = c("right", "left"), child = FALSE) {
  if(!isCount(degree, minimum = 0))
    stop2("`degree` must be a nonnegative integer: ", degree)
  if(!isCount(removal, minimum = 0))
    stop2("`removal` must be a nonnegative integer: ", removal)

  deg_right = deg_left = degree
  switch(match.arg(side),
         right = {deg_right <- deg_right + removal},
         left  = {deg_left <- deg_left + removal})

  # Chain on the left side
  x = linearPed(deg_left + 1)

  # Chain on the right side
  y = linearPed(deg_right + 1)
  y = relabel(y, old = 2:pedsize(y), new = pedsize(x) + 1:(pedsize(y)-1))

  # Merge
  z = mergePed(x, y)

  if (child) {
    parents = leaves(z)
    z = swapSex(z, parents[2], verbose = FALSE)
    z = addChildren(z, father = parents[1], mother = parents[2], nch = 1, verbose = FALSE)
  }
  z
}

#' @rdname ped_basic
#' @export
ancestralPed = function(g) {
  if(!isCount(g, minimum = 0))
    stop2("`g` must be a nonnegative integer: ", g)

  if(g == 0)
    return(singleton(1))

  N = 2^(g+1) - 1
  fathers = seq(from = 1, to = N-2, by = 2)
  mothers = fathers + 1

  ped(id = 1:N,
      fid = c(rep_len(0, 2^g), fathers),
      mid = c(rep_len(0, 2^g), mothers),
      sex = rep_len(1:2, N),
      reorder = FALSE, validate = FALSE, verbose = FALSE)
}

#' @rdname ped_basic
#' @export
selfingPed = function(s, sex = 1) {
  if(!isCount(s, minimum = 0))
    stop2("`s` must be a nonnegative integer: ", s)

  if(s == 0)
    return(singleton(1))

  ped(id = 1:(s+1), fid = 0:s, mid = 0:s, sex = c(rep(0, s), sex),
      reorder = FALSE, validate = FALSE, verbose = FALSE)
}

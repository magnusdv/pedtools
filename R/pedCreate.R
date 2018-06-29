#' Create simple pedigrees
#'
#' These are utility functions for creating some common pedigree structures as
#' `ped` objects. Use [swapSex()] to change the gender of
#' pedigree members.
#'
#' The call `cousinsPed(degree=n, removal=k)` creates a pedigree with two
#' n'th cousins, k times removed. By default, removals are added on the right
#' side. To override this, the parameter `degree2` can be used to indicate
#' explicitly the number of generations on the right side of the pedigree. When
#' `degree2` is given `removal` is ignored. (Similarly for
#' `halfCousinsPed`.)
#'
#' @param nch A positive integer indicating the number of offspring. If NULL, it
#'   is taken to be the `length(children)`
#' @param sex A numeric vector of length `nch` (recycled if shorter)
#'   encoding the genders of the children (0=unknown, 1=male, 2=female).
#' @param father The label of the father.
#' @param mother The label of the father.
#' @param children A character of length `nch`, with labels of the
#'   children.
#'
#' @param degree,degree2 Non-negative integers, indicating the degree of
#'   cousin-like relationships: 0=siblings, 1=first cousins; 2=second cousins,
#'   a.s.o. See Details and Examples.
#' @param removal Non-negative integers, indicating the removals of cousin-like
#'   relationships. See Details and Examples.
#' @param child A logical: Should an inbred child be added to the two cousins?
#'
#' @return A [ped()] object.
#'
#' @seealso [swapSex()], [removeIndividuals()],
#'   [addChildren()], [relabel()],
#'   [doubleCousins()],
#'
#' @examples
#'
#' # A nuclear family with 2 boys and 3 girls,
#' x = nuclearPed(5, sex=c(1,1,2,2,2))
#'
#' # Half sibs:
#' halfCousinsPed(degree=0)
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
#' @name pedCreate
NULL

#' @rdname pedCreate
#' @export
nuclearPed = function(nch, sex = 1, father = '1', mother = '2',
                      children = as.character(seq.int(3, length.out=nch))) {
  if(missing(nch))
    nch = length(children)
  assert_that(is.count(nch), length(children)==nch, length(father)==1,
              length(mother)==1, is.numeric(sex), length(sex) <= nch)

  if(length(sex) == 1) sex = rep(sex, nch)
  x = ped(id = 1:(2 + nch),
      fid = c(0, 0, rep(1, nch)),
      mid = c(0, 0, rep(2, nch)),
      sex = c(1, 2, sex))

  setLabels(x, c(father, mother, children))
}

#' @rdname pedCreate
#' @export
cousinsPed = function(degree, removal = 0, degree2 = NULL, child = FALSE) {
  # Creates a pedigree linking two n'th cousins k times removed, where
  # n=degree, k=removal. By default, removals are added on the right side,
  # i.e. degree2 = degree + removal.  If degree2 is non-NULL, the removal
  # parameter is ignored.
  assert_that(is_count0(degree), is_count0(removal), is.null(degree2) || is_count0(degree2))
  if (is.null(degree2))
      degree2 = degree + removal

  # Chain on the left side
  x = nuclearPed(1)
  for (i in seq_len(degree)) x = addSon(x, x$NIND, verbose = F)

  # Chain on the right side
  x = addChildren(x, father = 1, mother = 2, nch = 1, verbose = F)
  for (i in seq_len(degree2)) x = addSon(x, x$NIND, verbose = F)
  x = swapSex(x, x$NIND, verbose = T)

  if (child) {
      cous = leaves(x)
      x = addChildren(x, father = cous[1], mother = cous[2], nch = 1, verbose = F)
  }
  x
}

#' @rdname pedCreate
#' @export
halfCousinsPed = function(degree, removal = 0, degree2 = NULL, child = FALSE) {
  # Creates a pedigree linking two n'th half cousins k times removed, where
  # n=degree, k=removal.  By default, removals are added on the right side,
  # i.e. degree2 = degree + removal.  If degree2 is non-NULL, the removal
  # parameter is ignored.
  assert_that(is_count0(degree), is_count0(removal), is.null(degree2) || is_count0(degree2))
  if (is.null(degree2))
    degree2 = degree + removal

  # Chain on the left side
  x = nuclearPed(1)
  for (i in seq_len(degree)) x = addSon(x, x$NIND, verbose = F)

  # Chain on the right side
  x = addSon(x, 1, verbose = F)
  for (i in seq_len(degree2)) x = addSon(x, x$NIND, verbose = F)
  x = swapSex(x, x$NIND, verbose = F)

  if (child) {
    cous = leaves(x)
    x = addChildren(x, father = cous[1], mother = cous[2], nch = 1, verbose = F)
  }
  x
}


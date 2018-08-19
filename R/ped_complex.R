#' Complex pedigree structures
#'
#' These functions create certain classes of pedigrees that are not
#' straightforward to construct starting from the simple structures described in
#' [ped_basic].
#'
#' The function `doubleCousins` returns a pedigree linking two individuals who
#' are simultaneous paternal and maternal cousins. More precisely, they are:
#'
#' * paternal (full or half) cousins of type (`degree1`, `removal1`)
#'
#' * maternal (full or half) cousins of type (`degree2`, `removal2`).
#'
#' For convenience, a wrapper `doubleFirstCousins` is provided for the most
#' common case, double first cousins.
#'
#' `quadHalfFirstCousins` produces a pedigree with quadruple half first cousins.
#'
#' `fullSibMating` crosses full sibs continuously for the indicated number of
#' generations.
#'
#' `halfSibStack` produces a breeding scheme where the two individuals in the
#' final generation are simultaneously half siblings and half n'th cousins,
#' where `n=1,...,generations`.
#'
#'
#' @param degree1,degree2,removal1,removal2 Nonnegative integers.
#' @param half1,half2 Logicals, indicating if the fathers (resp mothers) should
#'   be full or half cousins.
#' @param child A logical: Should a child be added to the double cousins?
#' @param generations A positive integer indicating the number of crossings.
#'
#' @return A [`ped`] object.
#'
#' @seealso [ped_basic]
#'
#' @examples
#'
#' # Consecutive brother-sister matings.
#' fullSibMating(3)
#'
#' # Simultaneous half siblings and half first cousins
#' halfSibStack(2)
#'
#' # Double first cousins
#' x = doubleFirstCousins()
#' plot(x)
#'
#' # Quadruple half first cousins
#' QHFC = quadHalfFirstCousins()
#' # plot(QHFC) # Weird plotting behaviour for this pedigree.
#'
#' # A double half cousins pedigree with inbred founders,
#' # with the same identity coefficients as QHFC above.
#' # (Requires the `ribd` package.)
#' \dontrun{
#' x = doubleCousins(degree1 = 1, removal1 = 1, half1 = TRUE,
#'                   degree2 = 0, removal2 = 1, half2 = TRUE)
#' founder_inbreeding(x, 1) = 3-2*sqrt(2)
#' founder_inbreeding(x, 4) = .5 * sqrt(2)
#' j1 = ribd::ibd_identity(x, leaves(x))
#' j2 = ribd::ibd_identity(QHFC, leaves(QHFC))
#' stopifnot(identical(j1, j2))
#' }
#' @name ped_complex
NULL

#' @rdname ped_complex
#' @export
doubleCousins = function(degree1, degree2, removal1 = 0, removal2 = 0, half1 = FALSE, half2 = FALSE, child = FALSE) {
  if(!is_count(degree1, minimum=0)) stop2("`degree1` must be a nonnegative integer: ", degree1)
  if(!is_count(degree2, minimum=0)) stop2("`degree2` must be a nonnegative integer: ", degree2)
  if(!is_count(removal1, minimum=0)) stop2("`removal1` must be a nonnegative integer: ", removal1)
  if(!is_count(removal2, minimum=0)) stop2("`removal2` must be a nonnegative integer: ", removal2)

  # edge case: "paternal and maternal full 0th cousins". Interpreting this as full sibs for now. Review!
  if(degree1 + removal1 + half1 + degree2 + removal2 + half2 == 0) {
    x = nuclearPed(2, sex=1:2)
    if (child)
      x = addChildren(x, father = 3, mother = 4, nch = 1)
    return(x)
  }

  if(degree2*2+removal2 > degree1*2+removal1) {
    tmp = degree2; degree2 = degree1; degree1 = tmp
    tmp = removal2; removal2 = removal1; removal1 = tmp
  }

  # Paternal part
  if(half1)
    x1 = halfCousinsPed(degree1, removal = removal1)
  else
    x1 = cousinsPed(degree1, removal = removal1)

  offs = leaves(x1)

  # Maternal part
  if(degree2 + removal2 == 0) {
    if(half2) stop2("Full 0th cousins = full sibs, but this is incosistent with the other line")
    fath2 = father(x1, offs[2])
    moth1 = mother(x1, offs[1])
    x1 = removeIndividuals(x1, offs[2])

    x = addChildren(x1, father=fath2, mother=moth1, nch=1, sex=2, ids=offs[2])
  }
  else{
    if(half2)
      x2 = halfCousinsPed(degree2, removal = removal2)
    else
      x2 = cousinsPed(degree2, removal = removal2)
    x2 = swapSex(x2, labels(x2))
    x2 = swapSex(x2, leaves(x2))
    x2 = relabel(x2, seq(from = pedsize(x1)+1, length.out = pedsize(x2)))
    x2 = relabel(x2, old=leaves(x2), new=offs)

    x2 = relabel(x2, old=parents(x2, offs[1]), new=parents(x1, offs[1]))
    x2 = relabel(x2, old=parents(x2, offs[2]), new=parents(x1, offs[2]))
    x = mergePed(x1, x2)
  }

  if (child)
    x = addChildren(x, father = offs[1], mother = offs[2], nch = 1)
  x
}

#' @rdname ped_complex
#' @export
doubleFirstCousins = function()
    # Wrapper for the most common case of doubleCousins()
    doubleCousins(1, 1)

#' @rdname ped_complex
#' @export
quadHalfFirstCousins = function() {
    # Creates quad half fist cousins pedigree. NB: Does not draw well!
    ped(id=1:10, fid=c(0,0,0,0,1,3,1,3,5,7), mid=c(0,0,0,0,2,4,4,2,6,8),
        sex=rep(1:2, 5))
}

#' @rdname ped_complex
#' @export
fullSibMating = function(generations) {
    # Creates a pedigree resulting from repeated brother-sister matings.
    if(!is_count(generations)) stop2("`generations` must be a positive integer")
    x = nuclearPed(2, 1:2)
    for (i in seq_len(generations)[-1])
      x = addChildren(x, father = 2*i - 1, mother = 2*i, nch = 2, sex = 1:2, verbose = F)
    x
}

#' @rdname ped_complex
#' @export
halfSibStack = function(generations) {
    # Creates pedigree resulting from a breeding scheme where each generation adds two half
    # brothers and a female founder.  These become the parents of the half brothers in the next
    # layer.
    if(!is_count(generations)) stop2("`generations` must be a positive integer")
    x = ped(id = 1:5, fid = c(0, 0, 0, 1, 2), mid = c(0, 0, 0, 3, 3),
            sex = c(1, 1, 2, 1, 1), verbose = FALSE)
    for (g in seq_len(generations)[-1]) {
        m = 3 * g
        x = addChildren(x, father = m - 2, mother = m, nch = 1, verbose = F)
        x = addChildren(x, father = m - 1, mother = m, nch = 1, verbose = F)
    }
    x
}

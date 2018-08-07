#' Complex pedigree structures
#'
#' These functions create certain classes of pedigrees that are not
#' straightforward to construct starting from the simple structures described in [ped_basic].
#'
#' The function `doubleCousins` returns a pedigree linking two individuals whose
#' fathers are cousins of type (`degree1`, `removal1`), while the mothers are
#' cousins of type (`degree2`, `removal2`). For convenience, a wrapper
#' `doubleFirstCousins` is provided for the most common case, double first
#' cousins.
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
#' x = quadHalfFirstCousins()
#' # plot(x) # Weird plotting behaviour for this pedigree.
#'
#' @name ped_complex
NULL

#' @rdname ped_complex
#' @export
doubleCousins = function(degree1, degree2, removal1 = 0, removal2 = 0, half1 = FALSE, half2 = FALSE, child = FALSE) {
  if(!is_count(degree1, minimum=0)) stop2("`degree1` must be a nonnegative integer: ", degree1)
  if(!is_count(degree2, minimum=0)) stop2("`degree2` must be a nonnegative integer: ", degree2)
  if(!is_count(removal1, minimum=0)) stop2("`removal1` must be a nonnegative integer: ", removal1)
  if(!is_count(removal2, minimum=0)) stop2("`removal2` must be a nonnegative integer: ", removal2)

  if(degree2*2+removal2 > degree1*2+removal1) {
    tmp = degree2; degree2 = degree1; degree1 = tmp
    tmp = removal2; removal2 = removal1; removal1 = tmp
  }

  # Paternal part
  if(half1)
    x1 = halfCousinsPed(degree1, removal = removal1)
  else
    x1 = cousinsPed(degree1, removal = removal1)
  offs1 = leaves(x1)
  parents1 = parents(x1, offs1)
  moth1 = parents1[3]
  moth2 = parents1[4]

  if(degree2 == 0 && removal2 == 0) {
    stop2("This particular case is not implemented yet!")
    x1$MID[x1$MID == moth2] = moth1
    ped = ped[ -moth2, ] #?
    if (child)
      x = addChildren(x, father = offs1[1], mother = offs1[2], nch = 1)
    return(x)
  }

  # Maternal part
  if(half2)
    x2 = halfCousinsPed(degree2, removal = removal2)
  else
    x2 = cousinsPed(degree2, removal = removal2)
  offs2 = leaves(x2)
  parents2 = parents(x2, offs2)[c(3,4,1,2)]
  x2 = swapSex(x2, parents2)

  # Merging and adding children
  pat = setdiff(x1$LABELS, c(parents1, offs1)) # paternal line above parents
  mat = setdiff(x2$LABELS, c(parents2, offs2)) # maternal line above parents
  n = length(pat) + length(mat)

  x1 = relabel(x1, old=c(pat, parents1, offs1), new = c(1:length(pat), (1:6) + n))
  x2 = relabel(x2, old=c(mat, parents2, offs2), new = c(1:length(mat) + length(pat), (1:6) + n))

  x = mergePed(x1, x2)
  if (child)
    x = addChildren(x, father = offs1[1], mother = offs1[2], nch = 1)
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

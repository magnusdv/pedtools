#' Complex pedigree structures
#'
#' Functions for creating a selection of pedigrees that are awkward to construct
#' from scratch or with the simple structures described in [ped_basic].
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
#' `fullSibMating` crosses full sibs consecutively `n` times.
#'
#' `halfSibStack` produces a breeding scheme where the two individuals in the
#' final generation are simultaneous half k'th cousins, for each `k =
#' 0,...,n-1`.
#'
#' `halfSibTriangle` produces a triangular pedigree in which every pair of
#' parents are half siblings.
#'
#' @param degree1,degree2,removal1,removal2 Nonnegative integers.
#' @param half1,half2 Logicals, indicating if the fathers (resp. mothers) should
#'   be full or half cousins.
#' @param child A logical: Should a child be added to the double cousins?
#' @param n A positive integer indicating the number of crossings.
#'
#' @return A [`ped`] object.
#'
#' @seealso [ped_basic]
#'
#' @examples
#'
#' # Consecutive brother-sister matings.
#' x = fullSibMating(2)
#' # plot(x)
#'
#' # Simultaneous half siblings and half first cousins
#' x = halfSibStack(2)
#' # plot(x)
#'
#' # Double first cousins
#' x = doubleFirstCousins()
#' # plot(x)
#'
#' # Quadruple half first cousins
#' x = quadHalfFirstCousins()
#' # plot(x) # Weird plotting behaviour for this pedigree.
#'
#' # Triangular half-sib pattern
#' x = halfSibTriangle(4)
#' # plot(x)
#'
#' @name ped_complex
NULL

#' @rdname ped_complex
#' @export
doubleCousins = function(degree1, degree2, removal1 = 0, removal2 = 0,
                         half1 = FALSE, half2 = FALSE, child = FALSE) {
  if(!isCount(degree1, minimum = 0)) stop2("`degree1` must be a nonnegative integer: ", degree1)
  if(!isCount(degree2, minimum = 0)) stop2("`degree2` must be a nonnegative integer: ", degree2)
  if(!isCount(removal1, minimum = 0)) stop2("`removal1` must be a nonnegative integer: ", removal1)
  if(!isCount(removal2, minimum = 0)) stop2("`removal2` must be a nonnegative integer: ", removal2)

  # Ensure paternal path is longest (otherwise swap)
  if(2*degree1 + removal1 < 2*degree2 + removal2) {
    tmp = degree2; degree2 = degree1; degree1 = tmp
    tmp = removal2; removal2 = removal1; removal1 = tmp
    tmp = half2; half2 = half1; half1 = tmp
  }

  # Corner case: Paternal and maternal 0th cousins.
  if(degree1 + removal1 + degree2 + removal2 == 0) {
    if(half1 == half2) {
      x = nuclearPed(2, sex = 1:2)
      if (child)
        x = addChildren(x, father = 3, mother = 4, nch = 1)
      return(x)
    }
    else
      stop2("Paternal and maternal lines are incompatible.")
  }

  # Paternal part
  if(half1)
    x1 = halfCousinPed(degree1, removal = removal1)
  else
    x1 = cousinPed(degree1, removal = removal1)

  offs = leaves(x1)

  # Maternal part
  if(degree2 + removal2 == 0) {
    if(!half2)
      stop2("Paternal and maternal lines are incompatible.\n",
            "Maternal full 0th cousins => full siblings => fathers are identical.")
    fath2 = father(x1, offs[2])
    moth1 = mother(x1, offs[1])
    x1 = removeIndividuals(x1, offs[2],verbose = FALSE)

    x = addChildren(x1, father = fath2, mother = moth1, nch = 1, sex = 2, ids = offs[2])
  }
  else{
    if(half2)
      x2 = halfCousinPed(degree2, removal = removal2)
    else
      x2 = cousinPed(degree2, removal = removal2)

    # Change sex for everyone except the cousins
    x2 = swapSex(x2, .mysetdiff(labels(x2), leaves(x2)))
    x2 = relabel(x2, seq(from = pedsize(x1) + 1, length.out = pedsize(x2)))
    x2 = relabel(x2, old = leaves(x2), new = offs)

    x2 = relabel(x2, old = parents(x2, offs[1]), new = parents(x1, offs[1]))
    x2 = relabel(x2, old = parents(x2, offs[2]), new = parents(x1, offs[2]))
    x = mergePed(x1, x2)
  }

  if (child) {
    pars = leaves(x)
    x = swapSex(x, pars[2])
    x = addChildren(x, father = pars[1], mother = pars[2], nch = 1)
  }

  # Relabel according to plot order
  x = relabel(x, "asPlot")

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
    ped(id = 1:10,
        fid = c(0, 0, 0, 0, 1, 3, 1, 3, 5, 7),
        mid = c(0, 0, 0, 0, 2, 4, 4, 2, 6, 8),
        sex = rep(1:2, 5))
}

quadSecondCousins = function(type = c("cyclic", "exchange")) {
  type = match.arg(type)

  x = switch(type,
  cyclic = {
    pedmat = matrix(c(
      1,0,0,1,
      2,0,0,2,
      3,0,0,1,
      4,0,0,2,
      5,0,0,1,
      6,0,0,2,
      7,0,0,1,
      8,0,0,2,
      9,1,2,1,
      10,7,8,2,
      11,3,4,1,
      12,5,6,2,
      13,5,6,1,
      14,7,8,2,
      15,3,4,1,
      16,1,2,2,
      17,9,10,1,
      18,11,12,2,
      19,13,14,1,
      20,15,16,2,
      21,17,18,1,
      22,19,20,1), byrow = T, ncol = 4)

    ped(id  = pedmat[,1],
        fid = pedmat[,2],
        mid = pedmat[,3],
        sex = pedmat[,4],
        validate = FALSE, reorder = FALSE, verbose = FALSE) |>
      relabel()
  },
  exchange = {
    part1 = doubleFirstCousins() |>
      addSon(c(9,11), verbose = FALSE) |>
      addSon(c(10,13), verbose = FALSE)

    part2 = doubleFirstCousins() |>
      addSon(c(9,11), verbose = FALSE) |>
      addSon(c(10,13), verbose = FALSE) |>
      swapSex(9:10, verbose = FALSE) |>
      relabel(old = c(9,11, 10,13), new = c(11,9, 13,10))

    mergePed(part1, part2, by = c(9:14), relabel = TRUE)
  })
}

#' @rdname ped_complex
#' @export
fullSibMating = function(n) {
  if(!isCount(n, minimum = 0))
    stop2("`n` must be a nonnegative integer")

  id = seq_len(2*(n + 2))
  mothers = 2 * seq_len(n + 1)
  fathers = mothers - 1

  ped(id = id,
      fid = c(0, 0, rep(fathers, each = 2)),
      mid = c(0, 0, rep(mothers, each = 2)),
      sex = rep(1:2, times = n+2),
      validate = FALSE, reorder = FALSE, verbose = FALSE)
}

#' @rdname ped_complex
#' @export
halfSibStack = function(n) {
    # Creates pedigree resulting from a breeding scheme where each generation adds two half
    # brothers and a female founder.  These become the parents of the half brothers in the next
    # layer.
    if(!isCount(n)) stop2("`n` must be a positive integer")
    x = ped(id = 1:5, fid = c(0, 0, 0, 1, 2), mid = c(0, 0, 0, 3, 3),
            sex = c(1, 1, 2, 1, 1), verbose = FALSE)
    for (g in seq_len(n)[-1]) {
        m = 3 * g
        x = addChildren(x, father = m - 2, mother = m, nch = 1, verbose = FALSE)
        x = addChildren(x, father = m - 1, mother = m, nch = 1, verbose = FALSE)
    }
    x
}

#' @param g A positive integer; the number of generations.
#' @rdname ped_complex
#' @export
halfSibTriangle = function(g) {
  # A breeding scheme producing a triangular pedigree where every pair of parents are half sibs.
  if(!isCount(g))
    stop2("`g` must be a positive integer")
  if(g == 1)
    return(singleton(1))

  N = choose(g + 1, 2)
  id = seq_len(N)
  sex = unlist(lapply(g:1, function(a) rep_len(c(1L, 2L), length.out = a)))
  rw = rep(1:g, g:1) # row number
  par1 = pmax(0, id - g + rw - 2)
  par2 = pmax(0, id - g + rw - 1)
  fid = ifelse(sex == 1, par1, par2)
  mid = ifelse(sex == 1, par2, par1)

  ped(id = id, fid = fid, mid = mid, sex = sex, verbose = FALSE)
}

context("ped loops")

break_silent = function(...) breakLoops(..., verbose=F)
tie_silent = function(...) tieLoops(..., verbose=F)


test_that("loops are detected correctly in fullSibMating peds", {
  expect_false(fullSibMating(0)$UNBROKEN_LOOPS)
  expect_true(fullSibMating(1)$UNBROKEN_LOOPS)
  expect_true(fullSibMating(2)$UNBROKEN_LOOPS)
})

test_that("loops are gone after breaking in fullSibMating", {
  x = breakLoops(fullSibMating(1), verbose=F)
  expect_false(x$UNBROKEN_LOOPS)

  x = breakLoops(fullSibMating(2), verbose=F)
  expect_false(x$UNBROKEN_LOOPS)
})

test_that("breaking one loop at a time is same as all at once", {
  x = fullSibMating(2)
  x35 = break_silent(break_silent(x, 3), 5)
  x53 = break_silent(break_silent(x, 5), 3)
  expect_identical(as.data.frame(x35), as.data.frame(x53))
  expect_identical(tie_silent(x35), x)

  y35 = break_silent(x, c(3, 5))
  y53 = break_silent(x, c(5, 3))
  expect_identical(as.data.frame(y35), as.data.frame(y53))
  expect_identical(tie_silent(y35), x)
})

test_that("breakLoops gives same result when breakers are explicitly given", {
  x = fullSibMating(1)
  expect_identical(breakLoops(x, verbose=F), breakLoops(x, findLoopBreakers(x), verbose=F))

  x = fullSibMating(2)
  expect_identical(breakLoops(x, verbose=F), breakLoops(x, findLoopBreakers(x), verbose=F))
})

test_that("ped->matrix->ped works with broken loops (fullsib)", {
  x = breakLoops(fullSibMating(1), verbose=F)
  expect_identical(restorePed(as.matrix(x)), x)

  x = breakLoops(fullSibMating(2), verbose=F)
  expect_identical(restorePed(as.matrix(x)), x)
})

test_that("tieLoops restores broken loops (fullSib)", {
  x1 = fullSibMating(1)
  expect_identical(x1, tie_silent(break_silent(x1)))

  x2 = fullSibMating(2)
  expect_identical(x2, tie_silent(break_silent(x2)))

  x3 = reorderPed(x2, 8:1)
  expect_identical(x3, tie_silent(break_silent(x3)))

  x4 = relabel(x3, letters[1:8])
  expect_identical(x4, tie_silent(break_silent(x4)))
})

test_that("tieLoops restores broken IDENTITY ped", {
  skip("Identity examples - takes long time")
  ex = read.table (system.file ("example", "ex.pedigree", package = "identity"))
  x = ped(id=ex[,1], fid=ex[,2], mid=ex[,3], sex=rep(0,nrow(ex)), reorder=F)
  expect_identical(x, tie_silent(break_silent(x)))
})

test_that("loop breaking commutes with relabelling", {
  x = fullSibMating(2)
  x_r = relabel(x, old=1:8, new=letters[1:8])
  x_r_b = break_silent(x_r)

  x_b = break_silent(x)
  x_b_r = relabel(x_b, old=1:8, new=letters[1:8])

  skip("TODO: fix loop_breakers when relabelling...in relabel()")
  expect_identical(x_r_b, x_b_r)
})


# TODO: other looped peds


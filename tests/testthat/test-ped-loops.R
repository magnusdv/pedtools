
break_silent = function(...) breakLoops(..., allowFounder = TRUE, allowRepeated = TRUE, verbose=F)
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

  expect_identical(x_r_b, x_b_r)
})

test_that("it is possible to relabel a loop breaker copy individual", {
  x = break_silent(fullSibMating(1), 3)
  x_r = relabel(x, old="=3", new="copy of 3")

  lab = c(1:3, "copy of 3", 4:6)
  expect_identical(labels(x_r), lab)
})



# June 2026-----------------------------------------------------------------------------------------

test_that("founders can be loop breakers", {
  x = halfCousinPed(0, child = T)

  auto = findLoopBreakers(x, score = c(`2` = 1), allowFounder = TRUE)
  expect_identical(as.character(auto[, 1]), "2")

  plan = matrix(c("2", "4"), nrow = 1, dimnames = list(NULL, c("lb", "child")))

  y1 = break_silent(x, plan)
  y2 = break_silent(x, 2)
  expect_identical(y1,y2)
  expect_identical(tie_silent(y1), x)
})

test_that("loop breakers can be repeated", {
  x = halfSibTriangle(4)

  score = setNames(rep(-Inf, pedsize(x)), labels(x))
  score[c("3", "6")] = 1

  plan = findLoopBreakers(x, score = score, allowFounder = TRUE, allowRepeated = TRUE)
  expect_true(anyDuplicated.default(plan[, "lb"]) > 0)

  y1 = break_silent(x, plan)
  y2 = break_silent(x, plan[,1])
  expect_false(y1$UNBROKEN_LOOPS)
  expect_identical(y1,y2)
  expect_identical(tie_silent(y1), x)
})

test_that("repeated loop breakers are also checked across calls", {
  x = halfSibTriangle(4)

  z1 = x |> break_silent(c(6,6,3))
  z2 = x |> break_silent(6) |> break_silent(c(6,3)) |> reorderPed(neworder = z1$ID)
  expect_identical(z1, z2)
})

test_that("relabel() handles repeated loop breakers", {
  x = halfSibTriangle(4)
  y = break_silent(x, score = c("3" = 1, "6" = 1))

  z = relabel(y, c("6" = "A"))
  expect_identical(z$ID, sub("6", "A", y$ID))
  expect_identical(z$LOOP_BREAKERS, y$LOOP_BREAKERS)

  expect_identical(tie_silent(z) |> relabel(c(A = "6")), x)
})




test_that("getSex() works with and without labels", {
  x = nuclearPed(1)
  expect_equal(getSex(x, 1), 1)
  expect_equal(getSex(x, 1:3), c(1,2,1))

  y = relabel(x, c('fa', 'mo', 'ch'))
  expect_equal(getSex(y, 'mo'), 2)
  expect_equal(getSex(y, labels(y)), c(1,2,1))

})

test_that("swapSex() works in circular wedding loop", {
  x = addChildren(addSon(halfSibPed(), 3), 6, 1, 1)
  y = addChildren(addSon(halfSibPed(type = "mat"), 3), 1, 6, 1)
  expect_equal(swapSex(x, 1), y)
})

test_that("swapSex() works in ped list", {
  x = list(nuclearPed(), singleton("a"), singleton("b"), nuclearPed(fa=4, mo=5, ch=6))
  y = list(swapSex(nuclearPed(),1, verb=F), singleton("a", sex = 2),
           singleton("b"), nuclearPed(fa=4, mo=5, ch=6, sex=2))

  expect_equal(swapSex(x, c(1, "a", 6), verb=F), y)
})

test_that("is.pedList() is FALSE for empty list", {
  expect_false(is.pedList(list()))
})

test_that("pedsize works", {
  expect_equal(pedsize(singleton(1)), 1)
  expect_equal(pedsize(nuclearPed(1)), 3)

  x = fullSibMating(1)
  expect_equal(pedsize(x), 6)
  expect_equal(pedsize(breakLoops(x, verbose=F)), 7)
})

test_that("selfing is detected", {
  expect_false(hasSelfing(singleton(1)))
  expect_false(hasSelfing(nuclearPed(1)))

  x = addChildren(singleton(1, sex=0), 1,1,1)
  expect_true(hasSelfing(x))
})

test_that("common ancestors are detected", {
  labs = c("fa", "mo", "boy")
  x = relabel(nuclearPed(1), labs)

  ans = matrix(TRUE, ncol=3, nrow=3, dimnames=list(labs,labs))
  ans['fa','mo'] = ans['mo','fa'] = FALSE

  expect_identical(hasCommonAncestor(x), ans)
})

test_that("internalID gives empty output on empty intput", {
  y = nuclearPed(1)
  expect_equal(internalID(y, character()), integer(0))
  expect_equal(internalID(y, numeric()), integer(0))
})

test_that("internalID gives sensible error messages", {
  y = nuclearPed(1)
  expect_error(internalID(y, "foo"), "Unknown ID label: foo")
  expect_error(internalID(y, c("foo", "bar")), "Unknown ID label: foo, bar")
  expect_error(internalID(y, c("foo", 1:3, "bar")), "Unknown ID label: foo, bar")
  expect_error(internalID(y, ""), "Unknown ID label:")
  expect_error(internalID(y, 0), "Unknown ID label: 0")
})

test_that("internalID works in ped lists", {
  y = list(nuclearPed(1), singleton("a"))
  expect_equal(internalID(y, "a"), data.frame(id = "a", comp = 2L, int = 1L))
  expect_error(internalID(y, "foo"), "Unknown ID label: foo")
  expect_equal(internalID(y, "foo", errorIfUnknown = FALSE),
               data.frame(id = "foo", comp = NA_integer_, int = NA_integer_))
})

test_that("mergePed() works in half sib example", {
  x = nuclearPed(1)
  y = relabel(x, c(4,2,5))
  z = mergePed(x,y)
  zz = ped(id=1:5, fid=c(0,0,1,0,4), mid=c(0,0,2,0,2), sex=c(1,2,1,1,1))
  expect_identical(z,zz)
})

test_that("mergePed() works with `by`", {
  expect_identical(
    mergePed(nuclearPed(), nuclearPed(), by = 1:3, relabel = T),
    nuclearPed())

  expect_identical(
    mergePed(nuclearPed(), nuclearPed(), by = 1:3, relabel = T),
    nuclearPed())

  expect_identical(
    mergePed(nuclearPed(), nuclearPed(fa=1, mo=4,ch=5), by = 1, relabel = T),
    halfSibPed())

  expect_identical(
    mergePed(nuclearPed(fa = "A", mo = "B"), nuclearPed(), by = c(A=1), relabel = T),
    halfSibPed())
})

test_that("mergePed() catches errors", {

  expect_error(
    mergePed(nuclearPed(), nuclearPed(), by = c("4" = 2)),
    "Unknown ID label in pedigree 1")

  expect_error(
    mergePed(nuclearPed(), nuclearPed(), by = c("1" = 4)),
    "Unknown ID label in pedigree 2")

  expect_error(
    mergePed(nuclearPed(), nuclearPed(), by = c("1" = 2)),
    "Gender mismatch")

  expect_error(
    mergePed(nuclearPed(), nuclearPed(), by = 3),
    "Parent mismatch")
})

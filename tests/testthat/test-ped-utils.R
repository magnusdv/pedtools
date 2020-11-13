
test_that("getSex() works with and without labels", {
  x = nuclearPed(1)
  expect_equal(getSex(x, 1), 1)
  expect_equal(getSex(x, 1:3), c(1,2,1))

  y = relabel(x, c('fa', 'mo', 'ch'))
  expect_equal(getSex(y, 'mo'), 2)
  expect_equal(getSex(y, labels(y)), c(1,2,1))

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

test_that("mergePed() works in half sib example", {
  x = nuclearPed(1)
  y = relabel(x, c(4,2,5))
  z = mergePed(x,y)
  zz = ped(id=1:5, fid=c(0,0,1,0,4), mid=c(0,0,2,0,2), sex=c(1,2,1,1,1))
  expect_identical(z,zz)
})

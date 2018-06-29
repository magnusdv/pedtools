context("various utils")

test_that("pedSize of singleton is 1", {
  expect_equal(pedSize(singleton(1)), 1)
})

test_that("pedSize works", {
  expect_equal(pedSize(nuclearPed(1)), 3)
  expect_equal(pedSize(fullSibMating(2)), 6)
})

test_that("internalID gives empty output on empty intput", {
  y = nuclearPed(1)
  expect_equal(internalID(y, character()), integer(0))
  expect_equal(internalID(y, numeric()), integer(0))
})

test_that("internalID gives sensible error messages", {
  y = nuclearPed(1)
  expect_error(internalID(y, "foo"), "Unknown member of y: foo")
  expect_error(internalID(y, c("foo", "bar")), "Unknown members of y: foo, bar")
  expect_error(internalID(y, c("foo", 1:3, "bar")), "Unknown members of y: foo, bar")
  expect_error(internalID(y, ""), "Unknown member of y:")
  expect_error(internalID(y, 0), "Unknown member of y: 0")
})

test_that("mergePed() works in half sib example", {
  x = nuclearPed(1)
  y = relabel(x, c(4,2,5))
  z = mergePed(x,y)
  zz = ped(id=1:5, fid=c(0,0,1,0,4), mid=c(0,0,2,0,2), sex=c(1,2,1,1,1))
  expect_identical(z,zz)
})

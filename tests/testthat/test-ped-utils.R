context("various utils")

test_that("setting/getting founder inbreeding works properly", {
  x = addSon(nuclearPed(1), 3, id="boy", verbose=F)

  expect_equal(founder_inbreeding(x), c(0, 0, 0))
  expect_equal(founder_inbreeding(x, named=TRUE), c('1'=0, '2'=0, '4'=0))

  founder_inbreeding(x, '4') = 1
  expect_equal(founder_inbreeding(x, 4), 1)

  founder_inbreeding(x) = c('1' = 0.5)
  expect_equal(founder_inbreeding(x, named=TRUE), c('1'=0.5, '2'=0, '4'=1))
  expect_equal(founder_inbreeding(x, c(4,1)), c(1, 0.5))
})

test_that("founder inbreeding is preserved under modifications", {
  x = nuclearPed(1)
  founder_inbreeding(x, "1") = 1

  x1 = relabel(x, old=1, new="a")
  expect_equal(founder_inbreeding(x1, "a"), 1)

  x2 = reorderPed(x, 3:1)
  expect_equal(founder_inbreeding(x2, 1), 1)

  x3 = restore_ped(as.matrix(x))
  expect_equal(founder_inbreeding(x3, 1), 1)

  x4 = addSon(x, 1, verbose=F)
  expect_equal(founder_inbreeding(x4, 1), 1)

  x5 = addParents(x, 2, verbose=F)
  expect_equal(founder_inbreeding(x5, 1), 1)

  x6 = swapSex(x, 1, verbose=F)
  expect_equal(founder_inbreeding(x5, 1), 1)
})

test_that("founder inbreeding is removed by removeIndividuals()", {
  x = addSon(nuclearPed(1), 3, verbose=F)
  founder_inbreeding(x, 4) = 1
  expect_silent(y <- removeIndividuals(x, 4, verbose=F))
  expect_equal(founder_inbreeding(y), c(0,0))
})

test_that("setting/getting founder inbreeding catches errors", {
  x = addSon(nuclearPed(1), 3, id="boy", verbose=F)

  expect_error(founder_inbreeding("boy"), "Input is not a `ped` object")
  expect_error(founder_inbreeding(x, 5), "Unknown ID label")
  expect_error(founder_inbreeding(x, "boy"), "Pedigree member is not a founder")

  expect_error({founder_inbreeding(x, 5) = 0}, "Unknown ID label")
  expect_error({founder_inbreeding(x, x$LABELS) = 0}, "Pedigree member is not a founder: 3, boy")
  expect_error({founder_inbreeding(x, 1) = c(0,1)}, "Replacement vector must have length")
  expect_error({founder_inbreeding(x, c(1,1,2)) = c(0,0,0)}, "Duplicated ID label")
  expect_error({founder_inbreeding(x, 1) = "a"}, "Inbreeding coefficients must be numeric")
  expect_error({founder_inbreeding(x, 1:2) = c(-1,2)},
               "Inbreeding coefficients must be in the interval")

  founder_inbreeding(x, 1) = 1
  expect_error(addParents(x, 1, verbose=F),
               "Individual with nonzero founder inbreeding is no longer a founder")

})

test_that("getSex() works with and without labels", {
  x = nuclearPed(1)
  expect_equal(getSex(x, 1), 1)
  expect_equal(getSex(x, 1:3), c(1,2,1))

  y = relabel(x, c('fa', 'mo', 'ch'))
  expect_equal(getSex(y, 'mo'), 2)
  expect_equal(getSex(y, y$LABELS), c(1,2,1))

})

test_that("is.pedList() is FALSE for empty list", {
  expect_false(is.pedList(list()))
})

test_that("pedsize of singleton is 1", {
  expect_equal(pedsize(singleton(1)), 1)
})

test_that("pedsize works", {
  expect_equal(pedsize(nuclearPed(1)), 3)
  expect_equal(pedsize(fullSibMating(2)), 6)
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

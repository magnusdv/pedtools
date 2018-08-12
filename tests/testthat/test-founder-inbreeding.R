context("founder inbreeding")

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
  expect_error({founder_inbreeding(x, labels(x)) = 0}, "Pedigree member is not a founder: 3, boy")
  expect_error({founder_inbreeding(x, 1) = c(0,1)}, "Replacement vector must have length")
  expect_error({founder_inbreeding(x, c(1,1,2)) = c(0,0,0)}, "Duplicated ID label")
  expect_error({founder_inbreeding(x, 1) = "a"}, "Inbreeding coefficients must be numeric")
  expect_error({founder_inbreeding(x, 1:2) = c(-1,2)},
               "Inbreeding coefficients must be in the interval")

  founder_inbreeding(x, 1) = 1
  expect_error(addParents(x, 1, verbose=F),
               "Individual with nonzero founder inbreeding is no longer a founder")

})

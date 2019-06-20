context("founder inbreeding")

test_that("setting/getting founder inbreeding works properly", {
  x = addSon(nuclearPed(1), 3, id="boy", verbose=F)

  expect_equal(founderInbreeding(x), c(0, 0, 0))
  expect_equal(founderInbreeding(x, named=TRUE), c('1'=0, '2'=0, '4'=0))

  founderInbreeding(x, '4') = 1
  expect_equal(founderInbreeding(x, 4), 1)

  founderInbreeding(x) = c('1' = 0.5)
  expect_equal(founderInbreeding(x, named=TRUE), c('1'=0.5, '2'=0, '4'=1))
  expect_equal(founderInbreeding(x, c(4,1)), c(1, 0.5))
})

test_that("founder inbreeding is preserved under modifications", {
  x = nuclearPed(1)
  founderInbreeding(x, "1") = 1

  x1 = relabel(x, old=1, new="a")
  expect_equal(founderInbreeding(x1, "a"), 1)

  x2 = reorderPed(x, 3:1)
  expect_equal(founderInbreeding(x2, 1), 1)

  x3 = restorePed(as.matrix(x))
  expect_equal(founderInbreeding(x3, 1), 1)

  x4 = addSon(x, 1, verbose=F)
  expect_equal(founderInbreeding(x4, 1), 1)

  x5 = addParents(x, 2, verbose=F)
  expect_equal(founderInbreeding(x5, 1), 1)

  x6 = swapSex(x, 1, verbose=F)
  expect_equal(founderInbreeding(x5, 1), 1)
})

test_that("founder inbreeding is removed by removeIndividuals()", {
  x = addSon(nuclearPed(1), 3, verbose=F)
  founderInbreeding(x, 4) = 1
  expect_silent(y <- removeIndividuals(x, 4, verbose=F))
  expect_equal(founderInbreeding(y), c(0,0))
})

test_that("setting/getting founder inbreeding catches errors", {
  x = addSon(nuclearPed(1), 3, id="boy", verbose=F)

  expect_error(founderInbreeding("boy"), "Input is not a `ped` object")
  expect_error(founderInbreeding(x, 5), "Unknown ID label")
  expect_error(founderInbreeding(x, "boy"), "Pedigree member is not a founder")

  expect_error({founderInbreeding(x, 5) = 0}, "Unknown ID label")
  expect_error({founderInbreeding(x, labels(x)) = 0}, "Pedigree member is not a founder: 3, boy")
  expect_error({founderInbreeding(x, 1) = c(0,1)}, "Replacement vector must have length")
  expect_error({founderInbreeding(x, c(1,1,2)) = c(0,0,0)}, "Duplicated ID label")
  expect_error({founderInbreeding(x, 1) = "a"}, "Inbreeding coefficients must be numeric")
  expect_error({founderInbreeding(x, 1:2) = c(-1,2)},
               "Inbreeding coefficients must be in the interval")

  founderInbreeding(x, 1) = 1
  expect_message(addParents(x, 1, verbose=F), "Warning: Autosomal founder inbreeding lost.")

})

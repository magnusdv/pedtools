test_that("nuclearPed() direct labelling works", {
  expect_identical(nuclearPed(1), nuclearPed(children='3'))
  expect_identical(nuclearPed(1, sex=2), nuclearPed(children='3', sex=2))
  expect_identical(relabel(nuclearPed(1), letters[3:1]),
                   nuclearPed(fa='c', mo='b', children='a'))
  expect_identical(relabel(nuclearPed(1), old=3, new="foo"),
                   nuclearPed(children="foo"))
})

test_that("nuclearPed() catches errors", {
  expect_error(nuclearPed(), 'argument "nch" is missing')
  expect_error(nuclearPed(0), '`nch` must be a positive integer: 0')
  expect_error(nuclearPed(nch = 1:2), "`nch` must be a positive integer")
  expect_error(nuclearPed(nch = 2, child = 3), "`children` must have length `nch`")

  expect_error(nuclearPed(fa = 1:2, nch = 1), "`father` must have length 1")
  expect_error(nuclearPed(mo = 1:2, nch = 1), "`mother` must have length 1")

  expect_error(nuclearPed(nch = 1, sex='a'), "Illegal gender code: a")
  expect_error(nuclearPed(nch = 1, sex=1:2), "`sex` is longer than the number of individuals")
  expect_error(nuclearPed(nch = 1, sex=integer(0)), "`sex` cannot be empty")

  expect_error(nuclearPed(fa = 'a', child = 'a'), "Duplicated ID label: a")
  expect_error(nuclearPed(child = c("b", "b")), "Duplicated ID label: b")

  expect_error(nuclearPed(child = 1), "please specify a different label for the father")
  expect_error(nuclearPed(mo = 1, nch = 1), "please specify a different label for the father")
  expect_error(nuclearPed(child = 2), "please specify a different label for the mother")
  expect_error(nuclearPed(fa = 2, nch = 1), "please specify a different label for the mother")

})

test_that("halfSibPed() has expected ordering", {
  expect_equal(labels(halfSibPed()), as.character(1:5))
})

test_that("halfSibPed() recycles sex1 and sex2", {
  expect_equal(halfSibPed(2,3,sex1=0,sex2=2:1),
               halfSibPed(2,3,sex1=c(0,0), sex2=c(2,1,2)))
})

test_that("halfSibPed() catches errors", {
  expect_error(halfSibPed(-1), "`nch1` must be a positive integer")
  expect_error(halfSibPed(0), "`nch1` must be a positive integer")
  expect_error(halfSibPed(1, 0), "`nch2` must be a positive integer")

  expect_error(halfSibPed(sex1 = 'a'), "Illegal gender code: a")
  expect_error(halfSibPed(sex1 = 1:2), "`sex1` is longer than the number of individuals")
  expect_error(halfSibPed(sex1 = integer(0)), "`sex1` cannot be empty")
})

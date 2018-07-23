context("ped construction")

test_that("nuclearPed() works by giving labels", {
  expect_identical(nuclearPed(1), nuclearPed(children='3'))
  expect_identical(nuclearPed(1, sex=2), nuclearPed(children='3', sex=2))
  expect_identical(relabel(nuclearPed(1), letters[3:1]),
                   nuclearPed(fa='c', mo='b', children='a'))
  expect_identical(relabel(nuclearPed(1), old=3, new="foo"),
                   nuclearPed(children="foo"))
})

test_that("nuclearPed() gives sensible error messages", {
  expect_error(nuclearPed(), 'argument "nch" is missing')
  expect_error(nuclearPed(0), 'nch is not a count')
})

test_that("simple ped", {
  x1 = ped(id=1:3, fid=c(0,0,1), mid=c(0,0,2), sex=c(1,2,1), reorder=T)
  x2 = ped(id=1:3, fid=c(0,0,1), mid=c(0,0,2), sex=c(1,2,1), famid="")
  x3 = ped(id=1:3, fid=c(0,0,1), mid=c(0,0,2), sex=c(1,2,1), reorder=F)
  expect_is(x1, "ped")
  expect_identical(x1,x2)
  expect_identical(x1,x3)
})


test_that("random ped", {
  x = randomPed(3, 3)
  expect_is(x, "ped")
  expect_equal(pedSize(x), 6)
})

test_that("singleton creation works as expected", {
  x = singleton(1)
  expect_is(x, "ped")
  expect_is(x, "singleton")

  expect_error(singleton(), 'argument "id" is missing')
})


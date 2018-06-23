context("ped construction")

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
  expect_equal(x$NIND, 6)
})



test_that("generateLabs() works", {
  expect_equal(generateLabs("1"), "2")
  expect_equal(generateLabs("2"), "1")
  expect_equal(generateLabs("1", avoid = "2"), "3")
  expect_equal(generateLabs("2", avoid = "1"), "3")
  expect_equal(generateLabs("a"), "1")
  expect_equal(generateLabs(c("2", "a"), n = 2), c("1", "3"))
})

test_that("mergePed() works in half sib example", {
  x = nuclearPed(1)
  y = relabel(x, c(4,2,5))
  z = mergePed(x,y)
  zz = ped(id=1:5, fid=c(0,0,1,0,4), mid=c(0,0,2,0,2), sex=c(1,2,1,1,1))
  expect_identical(z,zz)
})

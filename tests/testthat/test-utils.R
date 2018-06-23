context("various utils")

test_that("nextNN works", {
  expect_equal(nextNN(c(1:5)), "NN_1")
  expect_equal(nextNN(character()), "NN_1")
  expect_equal(nextNN("NN1"), "NN_2")
  expect_equal(nextNN("NN.1"), "NN_2")
  expect_equal(nextNN("NN_1"), "NN_2")
  expect_equal(nextNN("NN-1"), "NN_2")
  expect_equal(nextNN(c(1,"NN_2", 3, "NN_1")), "NN_3")
})

test_that("mergePed() works in half sib example", {
  x = nuclearPed(1)
  y = relabel(x, c(4,2,5))
  z = mergePed(x,y)
  zz = ped(id=1:5, fid=c(0,0,1,0,4), mid=c(0,0,2,0,2), sex=c(1,2,1,1,1))
  expect_identical(z,zz)
})

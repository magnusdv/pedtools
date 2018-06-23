context("ped parts")

test_that("leaves() works with trivial labels", {
  expect_equal(leaves(nuclearPed(1)), "3")
  expect_equal(leaves(nuclearPed(8)), as.character(3:10))
})

test_that("leaves() works with labels", {
  x  = relabel(nuclearPed(1), c("f","m","c1"))
  expect_equal(leaves(x), 'c1')
  expect_equal(leaves(reorder(x, 3:1)), 'c1')
})

test_that("offspring are correct in nuclear", {
  x = nuclearPed(4)
  true_offs = as.character(3:6)
  expect_setequal(offspring(x, 1), true_offs)
  expect_setequal(offspring(x, 2), true_offs)
})

test_that("offspring are correct after reorder", {
  x = reorder(nuclearPed(4), 6:1)
  true_offs = as.character(3:6)
  expect_setequal(offspring(x, 1), true_offs)
  expect_setequal(offspring(x, 2), true_offs)
})

test_that("spouses are correct in nuclear", {
  x = nuclearPed(4)
  expect_equal(spouses(x,1), "2")
  expect_equal(spouses(x,2), "1")
  expect_equal(spouses(x,3), character(0))
  expect_equal(spouses(x,3, internal=T), numeric(0))
})

  test_that("spouses are correct after reorder", {
  x = reorder(nuclearPed(4), 6:1)
  expect_equal(spouses(x,1), "2")
  expect_equal(spouses(x,2), "1")
  expect_equal(spouses(x,3), character(0))
  expect_equal(spouses(x,3, internal=T), numeric(0))
})


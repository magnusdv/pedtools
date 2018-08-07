context("singletons")

test_that("singletons are constructed", {
  expect_error(singleton(1), NA)
  expect_error(singleton("A", sex=2, famid="FAM"), NA)
})

test_that("singleton modifications commute", {
  x = setLabels(swapSex(singleton(1), 1), 10)
  y = swapSex(relabel(singleton(1), 10), 10)
  expect_identical(x,y)

  z = relabel(singleton(1), old=1, new="A")
  expect_identical(z, singleton("A"))
})

test_that("addChildren works on singleton", {
  x = addChildren(singleton(1), father=1)
  expect_identical(x, nuclearPed(1))
})

test_that("addParents works on singleton", {
  x = addParents(singleton(3), 3, father=1, mother=2)
  x = parents_before_children(x)
  expect_identical(x, nuclearPed(1))
})

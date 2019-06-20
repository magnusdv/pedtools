context("singletons")

s = singleton(1)

test_that("singletons are constructed", {
  expect_error(s, NA)
  expect_error(singleton("A", sex=2, famid="FAM"), NA)
})

test_that("singleton modifications commute", {
  x = relabel(swapSex(s, 1), 10)
  y = swapSex(relabel(s, 10), 10)
  expect_identical(x,y)

  z = relabel(s, old=1, new="A")
  expect_identical(z, singleton("A"))
})

test_that("addChildren works on singleton", {
  x = addChildren(s, father=1)
  expect_identical(x, nuclearPed(1))
})

test_that("addParents works on singleton", {
  x = addParents(singleton(3), 3, father=1, mother=2)
  x = parentsBeforeChildren(x)
  expect_identical(x, nuclearPed(1))
})

test_that("singleton relatives are all empty", {
  expect_identical(father(s, 1), character(0))
  expect_identical(mother(s, 1), character(0))
  expect_identical(parents(s, 1), character(0))
  expect_identical(grandparents(s, 1), character(0))
  expect_identical(spouses(s, 1), character(0))
  expect_identical(children(s, 1), character(0))
  expect_identical(siblings(s, 1), character(0))
  expect_identical(cousins(s, 1), character(0))
  expect_identical(cousins(s, 1, half = T), character(0))
  expect_identical(cousins(s, 1, deg = 0, removal = 1), character(0))
  expect_identical(unrelated(s, 1), character(0))
  expect_identical(ancestors(s, 1), character(0))
  expect_identical(descendants(s, 1), character(0))
})

test_that("singleton subgroups are as expected", {
  expect_identical(leaves(s), "1")
  expect_identical(founders(s), "1")
  expect_identical(nonfounders(s), character(0))
  expect_identical(males(s), "1")
  expect_identical(females(s), character(0))
  expect_identical(typedMembers(s), character(0))
  expect_identical(untypedMembers(s), "1")
})

test_that("singleton with marker plots without error", {
  x = singleton(1)
  m = marker(x)
  pdf(NULL) # avoid Rplots.pdf being created
  expect_silent(plot(x, m))
})

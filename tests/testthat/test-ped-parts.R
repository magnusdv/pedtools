context("ped parts")

test_that("typedMembers() and untypedMembers() works", {
  x = nuclearPed(father="fa", mother="mo", child="ch")

  # labels
  expect_equal(typedMembers(x), character(0))
  expect_equal(untypedMembers(x), labels(x))

  # internal
  expect_equal(typedMembers(x, internal=T), numeric(0))
  expect_equal(untypedMembers(x, internal=T), 1:3)

  ### Add marker
  m = marker(x, 'mo'=1:2)
  y = setMarkers(x, m)

  # labels
  expect_equal(typedMembers(y), "mo")
  expect_equal(untypedMembers(y), c("fa", "ch"))

  # internal
  expect_equal(typedMembers(y, internal=T), 2)
  expect_equal(untypedMembers(y, internal=T), c(1,3))
})

test_that("leaves() works with trivial labels", {
  expect_equal(leaves(nuclearPed(1)), "3")
  expect_equal(leaves(nuclearPed(8)), as.character(3:10))
})

test_that("leaves() works with labels", {
  x  = relabel(nuclearPed(1), c("f","m","c1"))
  expect_equal(leaves(x), 'c1')
  expect_equal(leaves(reorderPed(x, 3:1)), 'c1')
})

test_that("internal = TRUE results in integer output", {
  x = singleton(1)
  expect_type(leaves(x, internal = T), "integer")
  expect_type(founders(x, internal = T), "integer")
  expect_type(nonfounders(x, internal = T), "integer")
  expect_type(males(x, internal = T), "integer")
  expect_type(females(x, internal = T), "integer")
  expect_type(typedMembers(x, internal = T), "integer")
  expect_type(untypedMembers(x, internal = T), "integer")
})


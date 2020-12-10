
test_that("father() and mother() works", {
  x = nuclearPed(father="fa", mother="mo", child="ch")
  expect_equal(father(x, "ch"), "fa")
  expect_equal(mother(x, "ch"), "mo")

  expect_equal(father(x, "fa"), character(0))
  expect_equal(mother(x, "fa"), character(0))

  expect_equal(father(x, 3, internal=T), 1)
  expect_equal(mother(x, 3, internal=T), 2)

  # TODO: these are not consistent with internal=F
  expect_equal(father(x, 1, internal=T), 0)
  expect_equal(mother(x, 1, internal=T), 0)
})


test_that("offspring are correct in nuclear", {
  x = nuclearPed(4)
  true_offs = as.character(3:6)
  expect_setequal(offspring(x, 1), true_offs)
  expect_setequal(offspring(x, 2), true_offs)
})

test_that("offspring are correct after reorder", {
  x = reorderPed(nuclearPed(4), 6:1)
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
  x = reorderPed(nuclearPed(4), 6:1)
  expect_equal(spouses(x,1), "2")
  expect_equal(spouses(x,2), "1")
  expect_equal(spouses(x,3), character(0))
  expect_equal(spouses(x,3, internal=T), numeric(0))
})

test_that("siblings are correctly detected", {
  x = halfSibPed(nch1 = 2, nch2 = 1)
  expect_equal(siblings(x, 1), character(0))
  expect_equal(siblings(x, 4), c("5", "6"))
  expect_equal(siblings(x, 4, half = T), "6")
  expect_equal(siblings(x, 4, half = F), "5")
})

test_that("unrelated ped members are correctly detected", {
  x = halfCousinPed(0, removal = 1, child = T)
  expect_equal(unrelated(x, 4), c("3", "6"))
  expect_equal(unrelated(x, 8), character(0))
})

test_that("internal = TRUE results in integer output", {
  x = singleton(1)
  expect_type(father(x, 1, internal = T), "integer")
  expect_type(mother(x, 1, internal = T), "integer")
  expect_type(parents(x, 1, internal = T), "integer")
  expect_type(children(x, 1, internal = T), "integer")
  expect_type(cousins(x, 1, internal = T), "integer")
  expect_type(ancestors(x, 1, internal = T), "integer")
  expect_type(descendants(x, 1, internal = T), "integer")
  expect_type(unrelated(x, 1, internal = T), "integer")
})

test_that("common ancestors are properly detected", {
  x = nuclearPed(2)
  expect_equal(commonAncestors(x, leaves(x)), as.character(1:2))
  expect_equal(commonAncestors(x, leaves(x), inclusive = T), as.character(1:2))
  expect_equal(commonAncestors(x, 2:3), character(0))
  expect_equal(commonAncestors(x, 2:3, inclusive = T), "2")
  expect_equal(commonAncestors(x, 1:4, inclusive = T), character(0))
})

test_that("common descendants are properly detected", {
  x = nuclearPed(2)
  expect_equal(commonDescendants(x, founders(x)), as.character(3:4))
  expect_equal(commonDescendants(x, founders(x), inclusive = T), as.character(3:4))
  expect_equal(commonDescendants(x, 2:3), character(0))
  expect_equal(commonDescendants(x, 2:3, inclusive = T), "3")
  expect_equal(commonDescendants(x, 1:4, inclusive = T), character(0))
})

test_that("common ancestors/descendants are detected in selfing ped", {
  x = selfingPed(3)
  expect_equal(commonAncestors(x, 1:4), character(0))
  expect_equal(commonAncestors(x, 1:4, inclusive = T), "1")
  expect_equal(commonAncestors(x, 2:4), "1")
  expect_equal(commonDescendants(x, 1:4), character(0))
  expect_equal(commonDescendants(x, 1:4, inclusive = T), "4")
  expect_equal(commonDescendants(x, 1:3), "4")
})

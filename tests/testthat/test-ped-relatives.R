
test_that("father() and mother() works", {
  x = nuclearPed(father="fa", mother="mo", children="ch")
  y = list(singleton(1), x, singleton(2))

  expect_equal(father(x, "ch"), "fa")
  expect_equal(mother(x, "ch"), "mo")
  expect_equal(father(y, "ch"), "fa")
  expect_equal(mother(y, "ch"), "mo")

  expect_equal(father(x, "fa"), character(0))
  expect_equal(mother(x, "fa"), character(0))
  expect_equal(father(y, "fa"), character(0))
  expect_equal(mother(y, "fa"), character(0))
  expect_equal(father(y, 1), character(0))
  expect_equal(mother(y, 1), character(0))

  expect_equal(father(x, 3, internal=T), 1)
  expect_equal(mother(x, 3, internal=T), 2)

  expect_equal(father(x, 1, internal=T), 0)
  expect_equal(mother(x, 1, internal=T), 0)

  expect_error(father(y, 1, internal=T),
               "Argument `internal` cannot be TRUE when `x` is disconnected")
  expect_error(mother(y, 1, internal=T),
               "Argument `internal` cannot be TRUE when `x` is disconnected")
})


test_that("children() are correct", {
  x = nuclearPed(4)
  y = list(singleton("NN1"), x, singleton("NN2"))

  true_offs = as.character(3:6)
  expect_setequal(children(x, 1), true_offs)
  expect_setequal(children(x, 2), true_offs)
  expect_setequal(children(y, 1), true_offs)
  expect_setequal(children(y, 2), true_offs)
})

test_that("children() are correct after reorder", {
  x = reorderPed(nuclearPed(4), 6:1)
  true_offs = as.character(3:6)
  expect_setequal(children(x, 1), true_offs)
  expect_setequal(children(x, 2), true_offs)
})


test_that("spouses() are correct", {
  x = nuclearPed(4)
  y = list(singleton("NN1"), x, singleton("NN2"))

  expect_equal(spouses(x,1), "2")
  expect_equal(spouses(x,2), "1")
  expect_equal(spouses(x,3), character(0))

  expect_equal(spouses(y,1), "2")
  expect_equal(spouses(y,2), "1")
  expect_equal(spouses(y,3), character(0))

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
  y = list(singleton("NN1"), x, singleton("NN2"))

  expect_equal(siblings(x, 1), character(0))
  expect_equal(siblings(x, 4), c("5", "6"))
  expect_equal(siblings(x, 4, half = T), "6")
  expect_equal(siblings(x, 4, half = F), "5")

  expect_equal(siblings(y, 1), character(0))
  expect_equal(siblings(y, 4), c("5", "6"))
  expect_equal(siblings(y, 4, half = T), "6")
  expect_equal(siblings(y, 4, half = F), "5")
})

test_that("unrelated ped members are correctly detected", {
  x = halfCousinPed(0, removal = 1, child = T)
  y = list(singleton("NN1"), x, singleton("NN2"))

  expect_equal(unrelated(x, 4), c("3", "6"))
  expect_equal(unrelated(x, 8), character(0))
  expect_equal(unrelated(y, 4), c("3", "6", "NN1", "NN2"))
  expect_equal(unrelated(y, 8), c("NN1", "NN2"))
})

test_that("internal = TRUE results in integer output", {
  x = singleton(1)
  expect_type(father(x, 1, internal = T), "integer")
  expect_type(mother(x, 1, internal = T), "integer")
  expect_type(parents(x, 1, internal = T), "integer")
  expect_type(children(x, 1, internal = T), "integer")
  # expect_type(cousins(x, 1, internal = T), "integer")
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

  y = list(singleton("NN1"), x, singleton("NN2"))
  expect_equal(commonAncestors(y, leaves(x)), as.character(1:2))
  expect_equal(commonAncestors(y, 2:3, inclusive = T), "2")
})

test_that("common descendants are properly detected", {
  x = nuclearPed(2)

  expect_equal(commonDescendants(x, founders(x)), as.character(3:4))
  expect_equal(commonDescendants(x, founders(x), inclusive = T), as.character(3:4))
  expect_equal(commonDescendants(x, 2:3), character(0))
  expect_equal(commonDescendants(x, 2:3, inclusive = T), "3")
  expect_equal(commonDescendants(x, 1:4, inclusive = T), character(0))

  y = list(singleton("NN1"), x, singleton("NN2"))
  expect_equal(commonDescendants(y, founders(x)), as.character(3:4))
  expect_equal(commonDescendants(y, 2:3, inclusive = T), "3")
})

test_that("common ancestors/descendants are detected in selfing ped", {
  x = selfingPed(3)
  expect_equal(commonAncestors(x, 1:4), character(0))
  expect_equal(commonAncestors(x, 1:4, inclusive = T), "1")
  expect_equal(commonAncestors(x, 2:4), "1")
  expect_equal(commonDescendants(x, 1:4), character(0))
  expect_equal(commonDescendants(x, 1:4, inclusive = T), "4")
  expect_equal(commonDescendants(x, 1:3), "4")

  y = list(singleton("NN1"), x, singleton("NN2"))
  expect_equal(commonAncestors(y, 1:4, inclusive = T), "1")
})

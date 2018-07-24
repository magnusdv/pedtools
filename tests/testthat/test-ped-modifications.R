context("ped modifications")

test_that("labels are correct after reordering", {
  x = setLabels(nuclearPed(1), letters[1:3])
  expect_identical(x$LABELS, letters[1:3])

  x = reorderPed(x, 3:1)
  expect_identical(x$LABELS, letters[3:1])
})

test_that("trivial reorder has no effect", {
  x = setLabels(nuclearPed(1), letters[1:3])
  y = reorderPed(x, 1:pedsize(x))
  expect_identical(x, y)
})

test_that("swapSex() works with trivial labels", {
  x = swapSex(nuclearPed(1, sex=1), 3)
  expect_equal(x$SEX[3], 2)
  expect_equal(getSex(x,3), 2)
})

test_that("swapSex() leaves parents unchanged", {
  x = nuclearPed(1)
  y = swapSex(x, 3)
  expect_equal(parents(x, 3), parents(y, 3))
})

test_that("swapSex() works with nonnumeric labels", {
  x = relabel(nuclearPed(1, sex=1), letters[1:3])
  y = swapSex(x, 'c')
  expect_equal(y$SEX[3], 2)
  expect_equal(getSex(y, 'c'), 2)
})

test_that("swapSex() works after reorder", {
  x = reorderPed(nuclearPed(1, sex=1), 3:1)
  y = swapSex(x, 3)
  expect_equal(y$SEX[1], 2)
  expect_equal(getSex(y, 3), 2)
})

test_that("swapSex() of all indivs", {
  x = swapSex(nuclearPed(1), 1:3)
  expect_equal(x$SEX, c(2,1,2))
  expect_equal(parents(x, 3), c("2", "1"))
})

test_that("addChildren works with num labels", {
  # start with male singleton
  m = singleton(4)
  expect_equal(addChildren(m, father=4, nch=2),
               nuclearPed(father=4, mother=5, children=6:7))
  expect_equal(addChildren(m, father=4, mother=1, nch=2, ids=2:3),
               nuclearPed(father=4, mother=1, children=2:3))

  # start with female singleton
  f = singleton(4, sex=2)

  f1 = addChildren(f, mother=4, nch=2)
  expect_equal(reorderPed(f1, c(2,1,3,4)),
               nuclearPed(mother=4, father=5, children=6:7))

  f2 = addChildren(f, mother=4, father=1, nch=2, ids=2:3)
  expect_equal(reorderPed(f2, c(2,1,3,4)),
               nuclearPed(mother=4, father=1, children=2:3))

  # inbreeding example
  x = nuclearPed(1, sex=2)
  y = addChildren(x, father=1, mother=3, nch=2, sex=2)
  expect_equal(spouses(y, 1), c("2", "3"))
  expect_equal(children(y, 1), c("3","4","5"))
  expect_equal(getSex(y, 4:5), c(2,2))
})

test_that("addChildren works with char labels", {
  m = singleton("fa")
  expect_equal(addChildren(m, father="fa", nch=2),
               nuclearPed(father="fa", mother="NN_1", children=c("NN_2","NN_3")))
  expect_equal(addChildren(m, father="fa", mother="mo", nch=2, ids=c("b1", "b2")),
               nuclearPed(father="fa", mother="mo", children=c("b1", "b2")))

  # start with female singleton
  f = singleton("mo", sex=2)

  f1 = addChildren(f, mother="mo", nch=2)
  expect_equal(reorderPed(f1, c(2,1,3,4)),
               nuclearPed(father="NN_1", mother="mo", children=c("NN_2","NN_3")))

  f2 = addChildren(f, mother="mo", father="fa", nch=2, ids=c("b1", "b2"))
  expect_equal(reorderPed(f2, c(2,1,3,4)),
               nuclearPed(father="fa", mother="mo", children=c("b1", "b2")))

  # inbreeding example
  x = nuclearPed(father="fa", mother="mo", children="da", sex=2)
  y = addChildren(x, father="fa", mother="da", nch=2, ids=c("g1", "g2"), sex=2)
  expect_equal(spouses(y, "fa"), c("mo", "da"))
  expect_equal(children(y, "fa"), c("da","g1","g2"))
  expect_equal(getSex(y, c("g1", "g2")), c(2,2))
})

test_that("addSon + addDaughter gives half sibs", {
  expect_equal(addDaughter(addSon(singleton(1), 1), 1),
               halfCousinsPed(0))
})

test_that("addChildren with nch=2 gives same result as twice with nch=1", {
  x = singleton(1)
  expect_equal(addChildren(addChildren(x, fa=1, mo=2), fa=1, mo=2),
               addChildren(x, fa=1, mo=2, nch=2))
})

test_that("adding and removing child restores original", {
  x = nuclearPed(1)
  y = addChildren(x, father=3, verbose=F)
  z = removeIndividuals(y, 5, verbose=F)
  expect_identical(x, z)

  x1 = fullSibMating(2)
  y1 = addChildren(x1, father=3, ids="99", verbose=F)
  z1 = removeIndividuals(y1, "99", verbose=F)
  expect_identical(x1, z1)

  x2 = setLabels(nuclearPed(1), c("F", "M", "C"))
  y2 = addChildren(x2, father="C", ids="baby", verbose=F)
  z2 = removeIndividuals(y2, "baby", verbose=F)
  expect_identical(x2, z2)

})

test_that("adding and removing child restores original - with markers", {
  x = nuclearPed(1)
  x = setMarkers(x, marker(x, '3'=1:2))
  y = addChildren(x, father=3, verbose=F)
  z = removeIndividuals(y, 5, verbose=F)
  expect_equal(x, z)
})

test_that("adding and removing parents restores original - with markers", {
  skip("`branch()` not implemented yet")
  x = nuclearPed(1)
  x = setMarkers(x, marker(x, '1'=1:2))
  y = addParents(x, id=1, verbose=F)
  z = branch(y, 1)
  expect_equal(x, z)
})

test_that("addParents() to nonfounder gives error", {
  x = nuclearPed(1)
  expect_error(addParents(x,3), "Individual 3 already has parents in the pedigree")

  y = relabel(x, letters[1:3])
  expect_error(addParents(y,'c'), "Individual c already has parents in the pedigree")
})

test_that("addParents() to multiple indivs gives error", {
  x = nuclearPed(1)
  expect_error(addParents(x, 1:2), "Only one individual at the time, please")
})

test_that("addParents() gives error if parent is impossible", {
  x = nuclearPed(1)
  expect_error(addParents(x, 1, father=3), "Assigned father is a descendant")
  expect_error(addParents(x, 1, mother=3), "Assigned mother is a descendant")
  y = addSon(x, 3, verbose=F)
  expect_error(addParents(y, 4, father=2), "Assigned father is female")
  expect_error(addParents(y, 4, mother=1), "Assigned mother is male")
})

test_that("addParents() gives message about new parents", {
  x = addSon(nuclearPed(1), 3, verbose=F)
  expect_message(addParents(x, 4, father=1), "Mother: Creating new individual with ID = 6")
  expect_message(addParents(x, 4, father=1, mother=123), "Mother: Creating new individual with ID = 123")
  expect_message(addParents(x, 4, father="1", mother="123"), "Mother: Creating new individual with ID = 123")

  expect_message(addParents(x, 4, mother=2), "Father: Creating new individual with ID = 6")
  expect_message(addParents(x, 4, father=123, mother=2), "Father: Creating new individual with ID = 123")
  expect_message(addParents(x, 4, father="123", mother="2"), "Father: Creating new individual with ID = 123")
})

test_that("addParents() creates new parents", {
  x = addSon(nuclearPed(1), 3, verbose=F)
  y1 = addParents(x, 4, verbose=F)
  y2 = addParents(x, 4, father=6, mother=7, verbose=F)
  expect_equal(y1$LABELS, as.character(1:7))
  expect_identical(y1, y2)
})

test_that("addParents() works with existing parents", {
  x = addSon(nuclearPed(1), 3, verbose=F)
  y = addParents(x, 4, father=1, mother=2, verbose=F)
  z = addChildren(nuclearPed(2, 1:2), 3, 4)
  expect_identical(y, z)
})


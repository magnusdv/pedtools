
test_that("labels are correct after reordering", {
  x = relabel(nuclearPed(1), letters[1:3])
  expect_identical(labels(x), letters[1:3])

  x = reorderPed(x, 3:1)
  expect_identical(labels(x), letters[3:1])
})

test_that("trivial reorder has no effect", {
  x = relabel(nuclearPed(1), letters[1:3])
  y = reorderPed(x, 1:pedsize(x))
  expect_identical(x, y)
})

test_that("setSex() reverses getSex()", {
  x = singleton(1)
  expect_identical(x, setSex(x, 1, getSex(x, 1)))
  expect_identical(x, setSex(x, sex = getSex(x, 1, named = T)))

  y = nuclearPed(fa="fa", child = "ch")
  expect_identical(y, setSex(y, labels(y), getSex(y)))
  expect_identical(y, setSex(y, sex = getSex(y, named = T)))

  z = singletons(c(4,1))
  expect_identical(z, setSex(z, sex = getSex(z, named = T)))
})

test_that("setSex() recycles sex", {
  x.male = singletons(1:2)
  x.female = singletons(1:2, sex = 2)
  expect_identical(setSex(x.male, ids = 1:2, sex = 2), x.female)

  y = nuclearPed(nch = 4)
  y2 = nuclearPed(nch = 4, sex = c(2,1,2,1))
  expect_identical(setSex(y, ids = 3:6, sex = 2:1), y2)
})

test_that("setSex() and swapSex() works with selectors", {
  x = nuclearPed(1) |> setSex(ids = leaves, sex = 0) |> swapSex(ids = males)
  y = nuclearPed(fa = 2, mo = 1, ch = 3, sex = 0) |> reorderPed()
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

test_that("swapSex() ingores indivs of unknown sex", {
  x = setSex(nuclearPed(1), ids = 3, sex = 0)
  y1 = swapSex(x, 3)
  expect_equal(getSex(y1, 3), 0)
  y2 = swapSex(x, 1:3)
  expect_equal(getSex(y2), c(2,1,0))
})

test_that("addChildren works with num labels", {
  # start with male singleton
  m = singleton(4) |> addChildren(father=4, nch=2)
  expect_equal(m, nuclearPed(father=4, mother=1, children=2:3))

  # start with female singleton
  f = singleton(4, sex=2) |> addChildren(mother=4, nch=2) |> reorderPed(c(1,4,2:3))
  expect_equal(f, nuclearPed(mother=4, father=1, children=2:3))

  # inbreeding example
  x = nuclearPed(1) |> addSon(2:3) |> addDaughter(2:3) |> addSon(2)
  expect_equal(spouses(x, 2), c("1", "3", "6"))
  expect_equal(children(x, 2), c("3","4","5","7"))
  expect_equal(getSex(x, leaves(x)), c(1,2,1))
})

test_that("addChildren works with char labels", {
  m = singleton("fa") |> addChildren(father="fa", nch=2)
  expect_equal(m, nuclearPed(father="fa", mother="1", children=2:3))

  # start with female singleton
  f = singleton("mo", sex=2) |> addChildren(mother="mo", nch=2) |> reorderPed(c(2,1,3,4))
  expect_equal(f, nuclearPed(father="1", mother="mo", children=2:3))
})

test_that("addChildren with nch=2 gives same result as twice with nch=1", {
  x = singleton(1)
  expect_equal(x |> addSon(1:2) |> addSon(1:2),
               addChildren(x, fa=1, mo=2, nch=2))

  expect_equal(x |> addSon(1:2, id = "B") |> addSon(1:2, id = "A"),
               addChildren(x, fa=1, mo=2, nch=2, ids = c("B", "A")))
})

test_that("adding and removing child restores original", {
  x = nuclearPed(1)
  y = x |> addSon(3, id = 5, verbose = F) |> removeIndividuals(5, verbose=F)
  expect_identical(x, y)

  x1 = fullSibMating(1)
  y1 = x1 |> addDaughter(3, id = "99", verbose=F) |> removeIndividuals("99", verbose=F)
  expect_identical(x1, y1)

  x2 = relabel(nuclearPed(1), c("F", "M", "C"))
  y2 = x2 |> addSon("C", id="baby", verbose=F) |> removeIndividuals("baby", verbose=F)
  expect_identical(x2, y2)

  # With marker
  xx = nuclearPed(1) |> addMarker('3' = "1/2")
  yy = xx |> addDaughter(3, id = 5, verbose=F) |> removeIndividuals(5, verbose=F)
  expect_equal(xx, yy)

})

test_that("removeIndividuals handles pedlists", {
  x = list(ancestralPed(2), singleton(8))
  y = removeIndividuals(x, 7) |> removeIndividuals(1:4, "ancestors") |>
    removeIndividuals(5:6)
  expect_identical(y, singleton(8))
})


test_that("addSon(), addDaughter(), addChild() creates children with correct sex", {
  x = singleton(1)
  expect_identical(x |> addSon(1, id = "A") |> getSex("A"), 1L)
  expect_identical(x |> addDaughter(1:2, id = 3) |> getSex(3), 2L)
  expect_identical(x |> addChild(c(1,"Mo"), id = "Ch", sex = 0) |> getSex("Ch"), 0L)
})

test_that("addSon() works with unordered parents", {
  x = nuclearPed(1)
  expect_identical(addSon(x, 1:2), addSon(x, 2:1))
  expect_identical(addSon(x, 3:4), addSon(x, 4:3))

  expect_identical(addDaughter(x, 1:2), addDaughter(x, 2:1))
  expect_identical(addDaughter(x, 3:4), addDaughter(x, 4:3))
})

test_that("addChildren() catches errors", {
  x = nuclearPed(1)
  expect_error(addChildren(x), "At least one parent must be an existing pedigree member")
  expect_error(addChildren(x, fa = 4), "At least one parent must be an existing pedigree member")
  expect_error(addChildren(x, mo = 4), "At least one parent must be an existing pedigree member")

  expect_error(addChildren(x, fa = 1:2), "More than one father indicated")
  expect_error(addChildren(x, mo = 1:2), "More than one mother indicated")

  expect_error(addChildren(x, 1, 2, nch = 0), "Argument `nch` must be a positive integer")
  expect_error(addChildren(x, 1, 2, nch = "a"), "Argument `nch` must be a positive integer")
  expect_error(addChildren(x, 1, 2, nch = 1:2), "Argument `nch` must be a positive integer")
  expect_error(addChildren(x, 1, 2, nch = list(1)), "Argument `nch` must be a positive integer")

  expect_error(addChildren(x, 1, 2, id = 3), "Individual already exist")
  expect_error(addChildren(x, 1, 2, id = 1:2), "Individual already exist")

  expect_error(addChildren(x, 1, 2, nch = 2, id = 1), "Length of `ids` must equal the number of children")
  expect_error(addChildren(x, 1, 2, nch = 1, id = 4:5), "Length of `ids` must equal the number of children")

  expect_error(addChildren(x, 1, 2, sex = -1), "Illegal value of `sex`")
  expect_error(addChildren(x, 1, 2, sex = NA), "Illegal value of `sex`")

  expect_error(addChildren(x, 1, 2, nch = 2, ids = c(4,4)), "Duplicated ID label")
})

test_that("addSon() and addDaughter() catches errors", {
  x = nuclearPed(1)
  expect_error(addSon(x, c(1,1)), "Duplicated parent")
  expect_error(addSon(x, c(1,3)), "Assigned mother is male")
  expect_error(addSon(x, 4:5), "At least one parent must be an existing pedigree member")

  expect_error(addDaughter(x, c(1,1)), "Duplicated parent")
  expect_error(addDaughter(x, c(1,3)), "Assigned mother is male")
  expect_error(addDaughter(x, 4:5), "At least one parent must be an existing pedigree member")
})

test_that("adding children across components", {
  x1 = singletons(1:3, sex = c(1,2,1)) |> addSon(1:2)
  x2 = list(nuclearPed(ch=4), singleton(3))
  expect_identical(x1, x2)

  y1 = singletons(c("a", "b", "d"), sex = c(1,2,2)) |>
    addSon(c("a", "b"), id = "c") |> addChildren("c", "d", id = "e", sex = 2)
  y2 = linearPed(2, sex = 1:2) |> relabel(letters[1:5])
  expect_identical(y1, y2)
})

test_that("modifaction chains give identical result", {
  x = singleton(3) |> addSon(3, id = "aa") |> addMarker(aa="1/1") |>
    addChild(c("aa", "bb"), id = "cc", sex = 0) |> setAlleleLabels(1, "A") |>
    relabel(c(cc = "c", bb = "b", aa = "a")) # |> plot(mark = 1)

  y = linearPed(2) |> setSex(5, sex = 0) |> addMarker(`3` = "A/-") |>
    relabel(c(3,1,"a","b","c")) |> setAlleles(ids = "a", marker = 1, alleles = "A")

  expect_identical(x,y)
})

test_that("adding and removing parents restores original - with markers", {
  x = nuclearPed(1) |> addMarker('1' = "1/2")
  y = addParents(x, id=1, verbose=F)
  z = branch(y, 1)
  expect_equal(x, z)
})

test_that("addParents() catches errors", {
  x = nuclearPed(1)
  expect_error(addParents(x,3),
               "Individual '3' already has parents in the pedigree")

  expect_error(x |> relabel(c("3" = "fa")) |> addParents("fa"),
               "Individual 'fa' already has parents in the pedigree")

  expect_error(addParents(x, 1:2), "Cannot add parents to multiple individuals")

  expect_error(addParents(x, 1, father=3), "Assigned father is a descendant")
  expect_error(addParents(x, 1, mother=3), "Assigned mother is a descendant")

  y = addSon(x, 3, verbose=F)
  expect_error(addParents(y, 4, father=2), "Assigned father is female")
  expect_error(addParents(y, 4, mother=1), "Assigned mother is male")
})

test_that("addParents() gives message about new parents", {
  x = addSon(nuclearPed(1), 3, verbose=F)
  expect_message(addParents(x, 4, father=1), "Creating new mother: 6")
  expect_message(addParents(x, 4, father=1, mother=123), "Creating new mother: 123")
  expect_message(addParents(x, 4, father="1", mother="123"), "Creating new mother: 123")

  expect_message(addParents(x, 4, mother=2), "Creating new father: 6")
  expect_message(addParents(x, 4, father=123, mother=2), "Creating new father: 123")
  expect_message(addParents(x, 4, father="123", mother="2"), "Creating new father: 123")
})

test_that("addParents() creates parents with correct labels", {
  x = nuclearPed(1)
  x1 = addParents(x, 1, verbose=F)
  x2 = addParents(x, 1, father=4, mother=5, verbose=F)
  expect_identical(x1, x2)

  y = nuclearPed(fa="fa", mo="mo", nch=1)
  y1 = addParents(y, "fa", verbose=F)
  y2 = addParents(y, "fa", father="1", mother="2", verbose=F)
  expect_identical(y1, y2)
})


test_that("addParents() works with existing parents", {
  x = addSon(nuclearPed(1), 3, verbose=F)
  y = addParents(x, 4, father=1, mother=2, verbose=F)
  z = addChildren(nuclearPed(2, 1:2), 3, 4)
  expect_identical(y, z)
})

test_that("addParents() adds parents before children", {
  x = addSon(nuclearPed(1), 3, verbose=F)
  x = addParents(x, 1, verbose=F)
  x = addParents(x, 4, father=1, verbose=F)
  expect_true(hasParentsBeforeChildren(x))
})

test_that("relabel() is strict", {
  x = nuclearPed(1)
  expect_error(relabel(x, old=c(3,3), new=4:5), "Duplicated entry in argument `old`: 3")
  expect_error(relabel(x, old=3, new=4:5), "Arguments `new` and `old` must have the same length")
  expect_error(relabel(x, old=3, new=2), "Duplicated ID label: 2")
})

test_that("relabelling is passed on to markers", {
  x = nuclearPed(1) |> addMarker()
  x = relabel(x, old=2, new="mother")
  expect_identical(attr(x$MARKERS[[1]], 'pedmembers'),
                   c('1','mother','3'))
})

test_that("relabel() works in pedlist", {
  x = list(nuclearPed(1), singleton(4))
  expect_identical(relabel(x, old=1, new="FA"),
                   list(nuclearPed(fa="FA", 1), singleton(4)))
  expect_identical(relabel(x, old=4, new="S"),
                   list(nuclearPed(1), singleton("S")))
  expect_identical(relabel(x, new=letters[1:4]),
                   list(nuclearPed(fa="a", mo = "b", ch = "c"), singleton("d")))
})

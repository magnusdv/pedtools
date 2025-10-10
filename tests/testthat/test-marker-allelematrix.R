
# Test case: singleton + trio, nonstandard labels.
df = data.frame(id = 4:1, fid = c(0,0,0,3), mid = c(0,0,0,2), sex = c(1,1,2,1),
                M1.1 = c(1,2,1,2), M1.2 = c(1,1,2,2), M2.1 = 3, M2.2 = 2)
x = as.ped(df, locus = list(list(alleles = 1:2), list(alleles = 2:3)))


test_that("`getAlleles()` returns a character matrix", {
  a = getAlleles(x, ids = "3", markers = "M1")
  expect_true(is.character(a) && is.matrix(a) && all(dim(a) == 1:2))

  b1 = getAlleles(x, ids = character(0), markers = "M1")
  expect_true(is.character(b1) && is.matrix(b1) && all(dim(b1) == c(0,2)))

  b2 = getAlleles(x, ids = "3", markers = character(0))
  expect_true(is.character(b2) && is.matrix(b2) && all(dim(b2) == c(1,0)))
})

test_that("`getAlleles()` and `setAlleles()` are inverses", {

  # All ids/markers
  y1 = setAlleles(x, alleles = getAlleles(x))
  expect_identical(x, y1)

  # Subset of ids
  ids = c(2,4)
  y2 = setAlleles(x, ids = ids, alleles = getAlleles(x, ids = ids))
  expect_identical(x, y2)

  # Subset of markers
  m = 2
  y3 = setAlleles(x, markers = m, alleles = getAlleles(x, markers = m))
  expect_identical(x, y3)

  # Subset of ids and markers
  ids = 3; m = 2
  y4 = setAlleles(x, ids = ids, markers = m, alleles = getAlleles(x, ids = ids, markers = m))
  expect_identical(x, y4)

  # Empty subset of ids/markers
  ids = character(); m = character(0)
  y5 = setAlleles(x, ids = ids, markers = m, alleles = getAlleles(x, ids = ids, markers = m))
  expect_identical(x, y5)
})

test_that("`setAlleles()` catches errors", {

  expect_error(setAlleles(x, alleles = 3), "Invalid allele for marker `M1`: 3")

  expect_error(setAlleles(x, ids = 4:5, alleles = 0), "Unknown ID label: 5")

  expect_error(setAlleles(x, markers = 3, alleles = 0), "Marker index out of range")

  expect_error(setAlleles(x, ids = 1, alleles = 1:3))
})

test_that("`removeGenotypes()` catches errors", {
  expect_error(removeGenotypes(x, ids = 4:5), "Unknown ID label: 5")
  expect_error(removeGenotypes(x, markers = 3), "Marker index out of range")
})

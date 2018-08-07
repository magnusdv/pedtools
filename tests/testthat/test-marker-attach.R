context("marker attachment")

test_that("setMarkers() creates empty markers with locus_annotations", {
  x = nuclearPed(1)
  ann = list(list(alleles=1:2, name="snp1"), list(alleles=c("a", "b")))
  x = setMarkers(x, allele_matrix = matrix(0, nrow=3, ncol=4), locus_annotations=ann)
  expect_equal(nMarkers(x), 2)
  expect_true(all(do.call(cbind, x$markerdata) == 0))
  expect_equal(name(x, 1), "snp1")
  expect_equal(alleles(x, 2), c("a", "b"))
})

test_that("addMarkers() adds empty markers with locus_annotations", {
  x = nuclearPed(1)
  x = setMarkers(x, marker(x, name="M0"))

  ann = list(list(alleles=1:2, name="snp1"), list(alleles=c("a", "b")))
  x = addMarkers(x, allele_matrix = matrix(0, nrow=3, ncol=4), locus_annotations=ann)
  expect_equal(nMarkers(x), 3)
  expect_equal(name(x, 2), "snp1")
  expect_equal(alleles(x, 3), c("a", "b"))
})


test_that("setMarkers() gives same result with and without allele_markers", {
  x = nuclearPed(fa="fa", mo="mo", children="boy")
  m = marker(x, mo='b', boy=c('b','c'), alleles=letters[1:3], name="snp")
  y = setMarkers(x, m)

  amatr = cbind(c(0,'b','b'),c(0,'b','c'))
  ann = list(alleles=letters[1:3], name="snp")
  z = setMarkers(x, allele_matrix = amatr, locus_annotations=ann)

  expect_identical(y,z)
})


test_that("addMarkers() gives same result with and without allele_markers", {
  x = nuclearPed(fa="fa", mo="mo", children="boy")
  x = setMarkers(x, marker(x, name="M0"))

  m = marker(x, mo='b', boy=c('b','c'), alleles=letters[1:3], name="snp")
  y = addMarkers(x, m)

  amatr = cbind(c(0,'b','b'),c(0,'b','c'))
  ann = list(alleles=letters[1:3], name="snp")
  z = addMarkers(x, allele_matrix = amatr, locus_annotations=ann)

  expect_identical(y,z)
})


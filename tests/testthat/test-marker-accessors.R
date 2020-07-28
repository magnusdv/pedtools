context("marker accessors")


test_that("simple marker getters work", {
  x = nuclearPed(1)
  m = marker(x, name="m1", chrom=1, posMb=1e7)
  x = setMarkers(x, m)

  expect_equal(name(m), "m1")
  expect_equal(name(x, 1), "m1")

  expect_equal(chrom(m), "1")
  expect_equal(chrom(x, markers=1), "1")
  expect_equal(chrom(x, markers="m1"), "1")

  expect_equal(posMb(m), 1e7)
  expect_equal(posMb(x, markers=1), 1e7)
  expect_equal(posMb(x, markers="m1"), 1e7)
})


test_that("alleles() accessor works", {
  x = nuclearPed(1)
  als = c("p","e","d")
  m1 = marker(x, alleles=1:3, name="m1")
  m2 = marker(x, alleles=als, name="m2")
  x = setMarkers(x, list(m1,m2))

  expect_equal(alleles(m1), as.character(1:3))
  expect_equal(alleles(x, marker=1), as.character(1:3))
  expect_equal(alleles(x, marker="m1"), as.character(1:3))

  expect_equal(alleles(m2), sort(als))
  expect_equal(alleles(x, marker=2), sort(als))
  expect_equal(alleles(x, marker="m2"), sort(als))
})

test_that("afreq() accessor works", {
  x = nuclearPed(1)
  afr = c(.2,.3,.5)
  m1 = marker(x, name="m1")
  m2 = marker(x, alleles=1:3, afreq=afr, name="m2")
  m3 = marker(x, alleles=3:1, afreq=afr, name="m3")
  x = setMarkers(x, list(m1,m2,m3))

  ans1 = c('1'=0.5, '2'=0.5)
  expect_equal(afreq(m1), ans1)
  expect_equal(afreq(x, marker=1), ans1)
  expect_equal(afreq(x, marker="m1"), ans1)

  names(afr) = 1:3
  expect_equal(afreq(m2), afr)
  expect_equal(afreq(x, marker=2), afr)

  afr_rev = rev(afr); names(afr_rev) = 1:3
  expect_equal(afreq(m3), afr_rev)
  expect_equal(afreq(x, marker=3), afr_rev)
})

test_that("afreq replacement works", {
  x = nuclearPed(1)
  m = marker(x, alleles=c("a", "b"), name="m1")
  x = setMarkers(x, list(m))

  afr = c(a=.1, b=.9)
  afreq(x, "m1") = afr
  expect_equal(afreq(x, 1), afr)

  afreq(x, 1) = rev(afr)
  expect_equal(afreq(x, "m1"), afr)

})

test_that("afreq replacement gives correct error messages", {
  x = nuclearPed(1)
  m = marker(x, alleles=c("a"), name="m1")
  x = setMarkers(x, list(m))

  expect_error({afreq(x, "m2") = c(a=1)}, "Unknown marker name: m2")
  expect_error({afreq(x, 2) = c(a=1)}, "Marker index out of range: 2")
  expect_error({afreq(x, 1:2) = c(a=1)}, "Frequency replacement can only be done for a single marker")
  expect_error({afreq(x, "m1") = 1}, "Frequency vector must be named")
  expect_error({afreq(x, "m1") = c(b=1)}, "Unknown allele: b")
  expect_error({afreq(x, "m1") = c(a=1)[0]}, "Alleles missing from frequency vector: a")
  expect_error({afreq(x, "m1") = c(a=0.1)}, "Frequencies must sum to 1")
})

test_that("genotype() works", {
  x = nuclearPed(children="boy") # labels are 1,2,boy
  m1 = marker(x, name="m1")
  m2 = marker(x, boy=1:2, name="m2")
  m3 = marker(x, '1'=17.2, name="m3") # homoz for STR allele
  x = setMarkers(x, list(m1,m2,m3))

  genoNA = c(NA_character_, NA_character_)
  expect_equal(genotype(m1, "boy"), genoNA)
  expect_equal(genotype(x, marker=1, id="boy"), genoNA)
  expect_equal(genotype(x, marker="m1", id="boy"), genoNA)

  genoHet = as.character(1:2)
  expect_equal(genotype(m2, id="boy"), genoHet)
  expect_equal(genotype(x, marker=2, id="boy"), genoHet)

  genoSTR = c("17.2", "17.2")
  expect_equal(genotype(m3, 1), genoSTR)
  expect_equal(genotype(m3, "1"), genoSTR)
  expect_equal(genotype(x, marker="m3", id=1), genoSTR)
})

test_that("genotype replacement works", {
  x = nuclearPed(father=101, mother=102, children="boy")
  m1 = marker(x, name="m1", alleles=1:2)
  m2 = marker(x, name="m2", alleles=c("a", "b"))
  x = setMarkers(x, list(m1, m2))

  genotype(x, 1, id=101) = 2
  genotype(x, "m1", "boy") = 2:1
  expect_equal(genotype(x, "m1", 101), c("2", "2"))
  expect_equal(genotype(x, 1, "boy"), c("2", "1"))

  genotype(x, 2, id=101) = 'b'
  genotype(x, "m2", "boy") = c('b','a')
  expect_equal(genotype(x, "m2", 101), c("b", "b"))
  expect_equal(genotype(x, 2, "boy"), c("b", "a"))
})

test_that("genotype replacement gives correct error messages", {
  x = nuclearPed(father=101, mother=102, children="boy")
  m1 = marker(x, name="m1", alleles=1:2)
  x = setMarkers(x, m1)

  expect_error({genotype(x, "m2", 101) = 3}, "Unknown marker name: m2")
  expect_error({genotype(x, 2, 101) = 3}, "Marker index out of range: 2")
  expect_error({genotype(x, 1:2, 101) = 3}, "Genotype replacement can only be done for a single marker")
  expect_error({genotype(x, "m1", 100) = 3}, "Unknown ID label: 100")
  expect_error({genotype(x, "m1", "girl") = 3}, "Unknown ID label: girl")
  expect_error({genotype(x, "m1", 101) = 3}, "Unknown allele for this marker: 3")
  expect_error({genotype(x, "m1", 101) = 1:3}, "Number of alleles must be 1 or 2")
})

test_that("genotype replacement works with partial genotypes", {
  x = nuclearPed(father=101, mother=102, children=1:2)
  m1 = marker(x, name="m1", alleles=c('a','b'))
  x = setMarkers(x, m1)

  genotype(x, "m1", id=101) = c("a", NA)
  genotype(x, "m1", id=102) = c("a", "")
  genotype(x, "m1", id=1) = c("b", 0)
  genotype(x, "m1", id=2) = c("b", "-")
  expect_equal(x$MARKERS[[1]][,1], c(1,1,2,2))
  expect_equal(x$MARKERS[[1]][,2], c(0,0,0,0))
  expect_equal(genotype(x, 1, 101), c("a", NA_character_))
  expect_equal(genotype(x, 1, 102), c("a", NA_character_))
  expect_equal(genotype(x, 1, 1), c("b", NA_character_))
  expect_equal(genotype(x, 1, 2), c("b", NA_character_))
})

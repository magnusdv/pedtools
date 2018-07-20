context("marker creation")

test_that("empty markers are created correctly", {
  x = nuclearPed(2)

  # Adding an empty SNP (all genotypes are missing):
  expect_is(marker(x, 0), "marker")
  expect_is(marker(x, '1'=0), "marker")
  expect_is(marker(x, '1'=1:2), "marker")

  # A SNP marker with alleles 'a' and 'b', for which
  # the father is heterozygous, the mother is homozygous for the
  # 'b' allele, first child has unknown genotype, while the last child is heterozygous.
  m1 = marker(x, '1'=c('a','b'), '2'='b', '4'=c('a','b'))

  # A rare SNP for which both children are heterozygous.
  # The 'alleles' argument can be skipped, but is recommended to ensure
  # correct order of the frequencies.
  m2 = marker(x, '3'=1:2, '4'=1:2, alleles=1:2, afreq=c(0.99, 0.01))


  # Similar shortcut for creating a marker for which all
  # pedigree members are homozygous for an allele (say 'b'):
  marker(x, 'b')
  # Alternative: m = marker(x, 'b'); addMarker(x, m)
})

test_that("alleles and freqs are sorted together", {
  x = nuclearPed(1)
  p = .8; q = 1-p

  # empty
  m1 = marker(x, alleles=2:1, afreq=c(q, p))
  expect_equal(attr(m1, 'alleles'), c('1','2'))
  expect_equal(attr(m1, 'afreq'), c(p, q))

  # the "2" allele observed
  m2 = marker(x, '1'=2, alleles=1:2, afreq=c(p, q))
  expect_equal(attr(m2, 'alleles'), c('1','2'))
  expect_equal(attr(m2, 'afreq'), c(p, q))

  # alleles not explicitly stated observed
  #m3 = marker(x, '1'=2:1, afreq=c(p, q))
  #expect_equal(attr(m3, 'alleles'), c('1','2'))
  #expect_equal(attr(m3, 'afreq'), c(p, q))
})

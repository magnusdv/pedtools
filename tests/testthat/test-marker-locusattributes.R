context("get/set global locus attributes")

# Test case 1: singleton + trivial markers
x = singleton(1)
x = addMarkers(x, marker(x, '1' = 1, name = "m1", chrom = "1"))
x = addMarkers(x, marker(x, '1' = 2, name = "m2", chrom = "X"))

# Test case 2: singleton + trio with two markers, nonstandard labels.
df = data.frame(id = 4:1,
                fid = c(0,0,0,3),
                mid = c(0,0,0,2),
                sex = c(1,1,2,1),
                A1 = c("1/1", "2/1", "1/2", "2/2"),
                A2 = "b/b")
attrs = list(list(name = "m1", alleles = 1:2),
             list(name = "m2", alleles = letters[1:2]))
y = as.ped(df, locusAttributes = attrs, sep = "/")

test_that("setLocusAttributes erases only when erase = T", {
  ll = list(alleles = 1:2)
  x1 = setLocusAttributes(x, markers = 1:2, locusA = ll, erase = T)
  x2 = setLocusAttributes(x, markers = 1:2, locusA = ll, erase = F)

  expect_true(all(is.na(name(x1, 1:2))))
  expect_true(all(is.na(chrom(x1, 1:2))))

  expect_identical(name(x2, 1:2), name(x, 1:2))
  expect_identical(chrom(x2, 1:2), chrom(x, 1:2))
})

test_that("setLocusAttributes catches errors", {
  sla = function(x = y, ...) setLocusAttributes(x, ...)

  expect_error(sla(marker = 2, loc = list(alleles = 1:2), erase = T),
               "Invalid allele for this marker")
  expect_error(sla(markers = 2, loc = list(alleles = 1:2), erase = F),
               "Invalid allele for marker `m2`")
  expect_error(sla(marker = 1, loc = list(alleles = 1, afreq = 2)),
               "Allele frequencies do not sum to 1")
  expect_error(sla(marker = 1, loc = list(afreq = c(.501, .501))),
               "Allele frequencies do not sum to 1")
  expect_error(sla(loc = list(name = "m1")),
               "When `locusAttributes` is a single list, then `markers` cannot be NULL")
  expect_error(sla(markers = 1:2, loc = list(list(name = 'm1'))),
               "List of locus attributes does not match the number of markers")
  expect_error(sla(setMarkers(x, NULL), loc = list()),
               "This function can only modify already attached markers")
})

test_that("setLocusAttributes matches on marker names", {
  sla = function(x = y, ...) setLocusAttributes(x, ...)

  locs = list(list(name = 'm1', alleles = 1:3),
              list(name = "m2", alleles = letters[1:3], afreq = c(.3, .3, .4)))

  expect_identical(sla(loc = locs[1], erase = T),
                   sla(markers = 1, loc = locs[1], erase = T))
  expect_identical(sla(loc = locs[1], erase = F),
                   sla(markers = 1, loc = locs[1], erase = F))

  expect_identical(sla(loc = locs[2:1]),
                   sla(markers = 1:2, loc = locs))
  expect_identical(sla(loc = locs[2:1], erase = T),
                   sla(markers = 1:2, loc = locs, erase = T))

  expect_identical(sla(loc = locs),
                   sla(markers = 2:1, loc = locs[2:1]))

})

test_that("set/get locusAttributes are inverses", {
  expect_identical(y, setLocusAttributes(y, loc = getLocusAttributes(y), erase = T))
  expect_identical(y, setLocusAttributes(y, markers = 2, loc = getLocusAttributes(y, markers = 2), erase = T))
})

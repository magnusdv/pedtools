
simpleMend = function(x, m) {
  x = setMarkers(x, m)
  mendelianCheck(x)
}

test_that("Autosomal mendelian errors are caught", {
  x = nuclearPed(father="fa", mother="mo", children=paste0("ch", 1:3))

  m1 = marker(x, 'fa'=1, 'mo'=1, 'ch1'=1:2)
  expect_message(simpleMend(x, m1), "Individual `ch1` incompatible with parents")

  m2 = marker(x, 'fa'=1, 'mo'=2, 'ch1'=1)
  expect_message(simpleMend(x, m2), "Individual `ch1` incompatible with parents")

  m3 = marker(x, 'fa'=1:2, 'mo'=1:2, 'ch1'=2:3)
  expect_message(simpleMend(x, m3), "Individual `ch1` incompatible with parents")

  m4 = marker(x, 'ch1'=1:2, 'ch2'=3:4, 'ch3'=5)
  expect_message(simpleMend(x, m4), "Sibship with parents `fa` and `mo` have too many alleles")

  m5 = marker(x, 'ch1'=1, 'ch2'=2:3, 'ch3'=3:4)
  expect_message(simpleMend(x, m5), "Sibship with parents `fa` and `mo` have too many alleles")

  m6 = marker(x, 'ch1'=1, 'ch2'=2, 'ch3'=3)
  expect_message(simpleMend(x, m6), "Sibship with parents `fa` and `mo` have too many alleles")
})

test_that("X-linked mendelian errors are caught", {
  x = nuclearPed(father="fa", mother="mo", children=paste0("ch", 1:3), sex=c(1,2,1))

  m1 = marker(x, 'fa'=2, 'mo'=1, 'ch1'=2, chrom=23)
  expect_message(simpleMend(x, m1), "Male `ch1` incompatible with mother")

  m2 = marker(x, 'fa'=1, 'mo'=2, 'ch2'=1, chrom=23)
  expect_message(simpleMend(x, m2), "Female `ch2` incompatible with parents")

  m3 = marker(x, 'fa'=1:2, chrom=23)
  expect_message(simpleMend(x, m3), "Male `fa` is heterozygous")

  skip("TODO: Mendelian check for X sibships")
  m4 = marker(x, 'ch1'=1, 'ch2'=3:4, 'ch3'=5, chrom=23)
  expect_message(simpleMend(x, m4), "Sibship with parents `fa` and `mo` have too many alleles")
})

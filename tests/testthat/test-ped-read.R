context("ped from data.frame or from file")

test_that("as.ped.data.frame identifies columns correctly", {
  A = singleton("A")
  dfA = data.frame(id="A", fid=0, mid=0, sex=1)
  expect_identical(as.ped(dfA), A)
  expect_identical(as.ped(rev(dfA)), A)

  B = singleton("B", sex=2, famid="1")
  dfB = data.frame(famid="1", id="B", fid=0, mid=0, sex=2)
  expect_identical(as.ped(dfB), B)
  expect_identical(as.ped(rev(dfB)), B)
})

test_that("data.frame with multiple singletons is converted to pedlist", {
  ans = list(S1 = singleton("A", sex=2, famid="S1"),
             S2 = singleton("B", sex=1, famid="S2"))
  df = data.frame(famid=c("S1", "S2"), id=c("A", "B"), fid=0, mid=0, sex=2:1)
  expect_identical(as.ped(df), ans)

  df_shuffled = rev(df)[2:1,]
  expect_identical(as.ped(df_shuffled), ans)
})

test_that("data.frame with multiple peds is converted to pedlist", {
  ans = list(S1 = nuclearPed(1),
             S2 = nuclearPed(1))
  famid(ans$S1) = "S1"
  famid(ans$S2) = "S2"

  df = data.frame(famid = rep(c("S1", "S2"), each=3),
                  id = rep(1:3,2), fid = rep(c(0,0,1), 2),
                  mid = rep(c(0,0,2), 2), sex = rep(c(1,2,1), 2))
  expect_identical(as.ped(df), ans)

  df_shuffled = rev(df)[6:1,]
  res = as.ped(df_shuffled)
  expect_identical(lapply(res, reorderPed), ans)
})

test_that("as.ped() converts data.frame with marker columns", {
  df = data.frame(id=c('fa','mo','boy'), fid=c(0,0,'fa'), mid=c(0,0,'mo'), sex=c(1,2,1),
                  c(0,0,1), c(0,0,2), c(0,0,2), c(0,0,2), fix.empty.names = F, stringsAsFactors = F)
  trio = nuclearPed(fa="fa", mo="mo", child="boy")
  x1 = as.ped(df)
  expect_equal(nMarkers(x1), 2)
  expect_identical(x1$markerdata[[1]], marker(trio, boy=1:2))
  expect_identical(x1$markerdata[[2]], marker(trio, boy=2))

  # force two alleles
  x2 = as.ped(df, locusAttributes = list(alleles=1:2))
  expect_identical(x2$markerdata[[1]], marker(trio, boy=1:2))
  expect_identical(x2$markerdata[[2]], marker(trio, boy=2, alleles=1:2))

  # Same, with markers in single columns:
  dfS = data.frame(id=c('fa','mo','boy'), fid=c(0,0,'fa'), mid=c(0,0,'mo'), sex=c(1,2,1),
                   c(NA,NA,"1/2"), c(NA,NA,"2/2"), fix.empty.names = F, stringsAsFactors = F)
  expect_identical(x1, as.ped(dfS, sep="/"))
  expect_identical(x2, as.ped(dfS, sep="/", locusAttributes = list(alleles=1:2)))
})

test_that("as.ped() converts data.frame to singletons with marker columns", {
  df = data.frame(famid=c("s1", "s2"), id=1, fid=0, mid=0, sex=1, m1=c(NA, "1/2"))
  pedlist = as.ped(df, sep="/")

  s1 = singleton(1, famid="s1")
  s1 = setMarkers(s1, marker(s1, name = "m1"))
  s2 = singleton(1, famid="s2")
  s2 = setMarkers(s2, marker(s2, '1'=1:2, name = "m1"))

  expect_identical(pedlist, list('s1'=s1, 's2'=s2))
})

test_that("as.ped() does not reorder (i.e. does not shuffle genotypes", {
  x = reorderPed(nuclearPed(1), 3:1)
  x = setMarkers(x, marker(x, '3' = 1:2, name = "M"))
  s = singleton("NN")
  s = setMarkers(s, marker(s, NN = 1:2, name = "M"))
  df = rbind(as.data.frame(x), as.data.frame(s))
  df = cbind(famid = c(1,1,1,2), df)

  y = as.ped(df, sep="/")
  famid(y[[1]]) = famid(y[[2]]) = ""

  expect_identical(list('1'=x, '2'=s), y)

})

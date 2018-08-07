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

test_that("clusterData works", {
  data("exps")

  expect_true(is.data.frame(exps))

  res <- clusterData(obj = exps,
                     clusterMethod = "kmeans",
                     clusterNum = 8)

  expect_true(is.list(res), "clusterData should return a list")

  expect_no_error(
    clusterData(
      obj = exps,
      clusterMethod = "kmeans",
      clusterNum = 8)
  )
})

test_that("visCluster works", {
  data("exps")

  expect_true(is.data.frame(exps))

  ck <- clusterData(obj = exps,
                    clusterMethod = "kmeans",
                    clusterNum = 8)

  expect_true(is.list(ck), "clusterData should return a list")

  p <- visCluster(object = ck,
                  plotType = "line")

  expect_s3_class(p, "ggplot")

  expect_no_error(
    visCluster(object = ck,
               plotType = "line")
  )
})

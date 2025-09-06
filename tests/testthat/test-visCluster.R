test_that("visCluster works", {
  data("exps")
  
  ck <- clusterData(
    obj = exps,
    cluster.method = "kmeans",
    cluster.num = 8)
  
  
  expect_no_error(
    visCluster(
      object = ck,
      plot.type = "line"
    )
  )
})

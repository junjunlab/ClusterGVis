test_that("clusterData works", {
  data("exps")
  
  expect_no_error(
    clusterData(
      obj = exps,
      cluster.method = "kmeans",
      cluster.num = 8)
  )
})

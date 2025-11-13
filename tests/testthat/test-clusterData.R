test_that("clusterData works", {
  data("exps")

  expect_true(is.data.frame(exps))

  expect_true(is.numeric(unlist(exps)))

  res <- clusterData(obj = exps,
                     clusterMethod = "kmeans",
                     clusterNum = 8)

  expect_true(is.list(res), "clusterData should return a list")

  # check output data content
  expect_true(all(c("wide.res","long.res","cluster.list","type","geneMode","geneType") %in% names(res)))

  expect_no_error(
    clusterData(
      obj = exps,
      clusterMethod = "kmeans",
      clusterNum = 8)
  )
})

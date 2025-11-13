test_that("visCluster works", {
  data("exps")

  expect_true(is.data.frame(exps))

  expect_true(is.numeric(unlist(exps)))

  ck <- clusterData(obj = exps,
                    clusterMethod = "kmeans",
                    clusterNum = 8)

  expect_true(is.list(ck), "clusterData should return a list")

  # check output data content
  expect_true(all(c("wide.res","long.res","cluster.list","type","geneMode","geneType") %in% names(ck)))

  p <- visCluster(object = ck,
                  plotType = "line")

  expect_true(ggplot2::is.ggplot(p))

  expect_no_error(
    visCluster(object = ck,
               plotType = "line")
  )

  p2 <- visCluster(object = ck,
                   plotType = "heatmap")

  expect_is(p2, "HeatmapList")

  expect_no_error(
    visCluster(object = ck,
               plotType = "heatmap")
  )
})

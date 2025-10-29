test_that("getClusters works", {
  data("exps")

  expect_true(is.data.frame(exps))

  result <- getClusters(obj = exps)

  expect_s3_class(result, "ggplot")

  expect_no_error(
    getClusters(obj = exps)
  )
})

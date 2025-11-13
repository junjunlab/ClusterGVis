test_that("getClusters works", {
  data("exps")

  expect_true(is.data.frame(exps))

  expect_true(is.numeric(unlist(exps)))

  result <- getClusters(obj = exps)

  expect_true(ggplot2::is.ggplot(result))

  expect_no_error(
    getClusters(obj = exps)
  )
})

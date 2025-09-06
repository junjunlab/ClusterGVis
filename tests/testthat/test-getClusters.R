test_that("getClusters works", {
  data("exps")
  expect_no_error(
    getClusters(obj = exps)
  )
})

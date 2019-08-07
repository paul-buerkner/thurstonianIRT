context("Tests for TIRT models fitted with Mplus")

skip("Mplus is not open source software")

test_that("mplus code for bernoulli responses works", {
  set.seed(1234)
  lambdas <- c(runif(6, 0.5, 1), runif(6, -1, -0.5))
  sdata <- sim_TIRT_data(
    npersons = 100,
    ntraits = 3,
    nblocks_per_trait = 4,
    gamma = 0,
    lambda = lambdas,
    Phi = diag(3),
    family = "bernoulli"
  )
  fit <- suppressWarnings(fit_TIRT_mplus(sdata))
  expect_is(fit, "TIRTfit")
  pr <- suppressWarnings(predict(fit))
  pr_names <- c("id", "trait", "estimate")
  expect_equal(names(pr), pr_names)
  expect_equal(length(unique(pr$id)), 100)
})

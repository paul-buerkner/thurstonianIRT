context("Tests for TIRT models fitted with lavaan")

test_that("lavaan code for bernoulli responses works", {
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
  fit <- suppressWarnings(fit_TIRT_lavaan(sdata))
  expect_is(fit, "TIRTfit")
  pr <- suppressWarnings(predict(fit))
  pr_names <- c("id", "trait", "estimate")
  expect_equal(names(pr), pr_names)
  expect_equal(length(unique(pr$id)), 100)

  # test predictions for new data
  new_sdata <- dplyr::filter(sdata, person %in% 1:5)
  pr_new <- suppressWarnings(predict(fit, new_sdata))
  expect_equal(names(pr_new), pr_names)
  expect_equal(length(unique(pr_new$id)), 5)
})

test_that("lavaan code for gaussian responses works", {
  set.seed(12345)
  lambdas <- c(runif(6, 0.5, 1), runif(6, -1, -0.5))
  sdata <- sim_TIRT_data(
    npersons = 100,
    ntraits = 3,
    nblocks_per_trait = 4,
    gamma = 0,
    lambda = lambdas,
    Phi = diag(3),
    family = "gaussian"
  )
  fit <- suppressWarnings(fit_TIRT_lavaan(sdata))
  expect_is(fit, "TIRTfit")
  pr <- suppressWarnings(predict(fit))
  pr_names <- c("id", "trait", "estimate")
  expect_equal(names(pr), pr_names)
  expect_equal(length(unique(pr$id)), 100)

  # test predictions for new data
  new_sdata <- dplyr::filter(sdata, person %in% 1:5)
  pr_new <- suppressWarnings(predict(fit, new_sdata))
  expect_equal(names(pr_new), pr_names)
  expect_equal(length(unique(pr_new$id)), 5)
})

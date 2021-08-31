context("Tests for TIRT models fitted with Stan")

test_that("Stan code for bernoulli responses works", {
  set.seed(1234)
  lambdas <- c(runif(6, 0.5, 1), runif(6, -1, -0.5))
  sdata <- sim_TIRT_data(
    npersons = 10,
    ntraits = 3,
    nblocks_per_trait = 4,
    gamma = 0,
    lambda = lambdas,
    Phi = diag(3),
    family = "bernoulli"
  )
  fit <- suppressWarnings(fit_TIRT_stan(sdata, chains = 1, iter = 500))
  expect_is(fit, "TIRTfit")
  pr <- predict(fit)
  pr_names <- c("id", "trait", "estimate", "se", "lower_ci", "upper_ci")
  expect_equal(names(pr), pr_names)
  expect_equal(length(unique(pr$id)), 10)

  # test predictions for new data
  new_sdata <- dplyr::filter(sdata, person %in% 1:5)
  pr_new <- predict(fit, new_sdata, chains = 1, iter = 500)
  expect_equal(names(pr_new), pr_names)
  expect_equal(length(unique(pr_new$id)), 5)
})

test_that("Stan code for ordinal responses works", {
  set.seed(1234)
  ncat <- 4
  gamma <- matrix(
    seq(-2, 2, length.out = max(ncat) - 1),
    nrow = 12,
    ncol = max(ncat) - 1,
    byrow = TRUE
  )
  lambdas <- c(runif(6, 0.5, 1), runif(6, -1, -0.5))
  sdata <- sim_TIRT_data(
    npersons = 10,
    ntraits = 3,
    nblocks_per_trait = 4,
    gamma = gamma,
    lambda = lambdas,
    Phi = diag(3),
    family = "cumulative"
  )
  fit <- suppressWarnings(fit_TIRT_stan(sdata, chains = 1, iter = 500))
  expect_is(fit, "TIRTfit")
  pr <- predict(fit)
  pr_names <- c("id", "trait", "estimate", "se", "lower_ci", "upper_ci")
  expect_equal(names(pr), pr_names)
  expect_equal(length(unique(pr$id)), 10)

  # test predictions for new data
  new_sdata <- dplyr::filter(sdata, person %in% 1:5)
  pr_new <- predict(fit, new_sdata, chains = 1, iter = 500)
  expect_equal(names(pr_new), pr_names)
  expect_equal(length(unique(pr_new$id)), 5)
})

test_that("Stan code for gaussian responses works", {
  set.seed(1234)
  lambdas <- c(runif(6, 0.5, 1), runif(6, -1, -0.5))
  sdata <- sim_TIRT_data(
    npersons = 10,
    ntraits = 3,
    nblocks_per_trait = 4,
    gamma = 0,
    lambda = lambdas,
    Phi = diag(3),
    family = "gaussian"
  )
  fit <- suppressWarnings(fit_TIRT_stan(sdata, chains = 1, iter = 500))
  expect_is(fit, "TIRTfit")
  pr <- predict(fit)
  pr_names <- c("id", "trait", "estimate", "se", "lower_ci", "upper_ci")
  expect_equal(names(pr), pr_names)
  expect_equal(length(unique(pr$id)), 10)

  # test predictions for new data
  new_sdata <- dplyr::filter(sdata, person %in% 1:5)
  pr_new <- predict(fit, new_sdata, chains = 1, iter = 500)
  expect_equal(names(pr_new), pr_names)
  expect_equal(length(unique(pr_new$id)), 5)
})

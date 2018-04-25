context("Tests for TIRT models fitted with lavaan")

test_that("lavaan code runs without errors", {
  lambdas <- c(runif(6, 0.5, 1), runif(6, -1, -0.5))
  sdata <- sim_thurstonian_data(
    npersons = 100,
    ntraits = 3,
    nblocks_per_trait = 4,
    gamma = 0,
    lambda = lambdas,
    Phi = diag(3)
  )
  lfit <- suppressWarnings(fit_TIRT_lavaan(sdata))
  expect_is(lfit, "TIRTfit")
})




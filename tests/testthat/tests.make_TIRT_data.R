context("Tests for creating T-IRT compatible data sets")

test_that("make_TIRT_data works correctly", {
  data("triplets")
  # define the blocks of items
  blocks <-
    set_block(c("i1", "i2", "i3"), traits = c("t1", "t2", "t3"),
              signs = c(1, 1, 1)) +
    set_block(c("i4", "i5", "i6"), traits = c("t1", "t2", "t3"),
              signs = c(-1, 1, 1)) +
    set_block(c("i7", "i8", "i9"), traits = c("t1", "t2", "t3"),
              signs = c(1, 1, -1)) +
    set_block(c("i10", "i11", "i12"), traits = c("t1", "t2", "t3"),
              signs = c(1, -1, 1))
  # generate the data to be understood by 'thurstonianIRT'
  triplets_long <- make_TIRT_data(
    data = triplets, blocks = blocks, direction = "larger",
    format = "pairwise", family = "bernoulli", range = c(0, 1)
  )
  expect_equal(NROW(triplets_long), 1200)
  expect_true(all(triplets_long$response %in% 0:1))
})

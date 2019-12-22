# these functions are likely redundant to make_TIRT_data(format = "ranks")

# compare two values; NA if equal
greater <- function(x) {
  stopifnot(length(x) == 2L)
  if (x[1] > x[2]) {
    out <- 1L
  } else if (x[1] < x[2]) {
    out <- 0L
  } else {
    out <- NA_integer_
  }
  out
}

#' Transform ranks into binary comparisons
#'
#' Transform data of ranks into data of binary item comparisons.
#' The output can be used in function \code{\link{make_TIRT_data}}.
#'
#' @param data A \code{data.frame} containing all ranks for every block.
#' @param nitems_per_block An integer indicating the number of items per block.
#' All items within the same block should be adjacent in \code{data}.
#'
#' @examples
#' # define ranking data with two rows (participants)
#' # and two blocks with three items each
#' ranks <- data.frame(
#'   i1 = c(1,2), i2 = c(2,3), i3 = c(3,1),
#'   i4 = c(3,1), i5 = c(1,3), i6 = c(2,2)
#' )
#' ranks
#'
#' # generate data of binary comparisons
#' transform_ranks(ranks, nitems_per_block = 3)
#'
#' @seealso \code{\link{make_TIRT_data}}
#' @noRd
transform_ranks <- function(data, nitems_per_block) {
  data <- as.data.frame(data)
  nblocks <- ncol(data) / nitems_per_block
  if (nitems_per_block %% 1 != 0 | nitems_per_block <= 1){
    stop("Number of items per block must be an integer greater than 1.")
  }
  if (nblocks %% 1 != 0){
    stop("Number of columns in 'data' must be dividable by 'nitems_per_block'.")
  }
  splitted <- split.default(
    data, rep(seq_len(nblocks), rep(nitems_per_block, nblocks))
  )
  splitted_dich <- lapply(splitted, transform_ranks_per_block)
  do.call(cbind, splitted_dich)
}

# create pairs of comparisons for the block
compare_item_pairs <- function(x) {
  combn(unlist(x), 2, greater)
}

# convert block of ranks to block of paired comparisons with appropriate column names
transform_ranks_per_block <- function(ranks) {
  out <- t(apply(ranks, 1, compare_item_pairs))
  colnames(out) <- combn(names(ranks), 2, collapse)
  out
}


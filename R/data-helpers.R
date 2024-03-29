#' Prepare data for Thurstonian IRT models
#'
#' @param data An object of class \code{data.frame} containing data of all
#'   variables used in the model.
#' @param blocks Object of class \code{TIRTblocks} generated by
#'   \code{\link{set_block}} indicating which items belong to which block, trait
#'   and more. Ignored if data already contains information on the blocks.
#' @param direction Indicates if \code{"larger"} (the default) or
#'   \code{"smaller"} input values are considered as indicating the favored
#'   answer.
#' @param format Format of the item responses. Either \code{"ranks"} for
#'   responses in ranked format or \code{"pairwise"} for responses in pairwise
#'   comparison format. If \code{"ranks"}, each item must have its own
#'   column in the data frame which contains its ranks within the block.
#'   If \code{"pairwise"}, each existing item combination must have its
#'   own column named after the combination of the two compared items.
#' @param family Name of assumed the response distribution. Either
#'   \code{"bernoulli"}, \code{"cumulative"}, or \code{"gaussian"}.
#' @param partial A flag to indicate whether partial comparisons are allowed
#'   for responses stored in the \code{"ranks"} format.
#' @param range Numeric vector of length two giving the range of the
#'   responses when using the \code{"pairwise"} format. Defaults
#'   to \code{c(0, 1)} for use with dichotomous responses.
#'
#' @return A \code{data.frame} in a specific format and with attributes ready
#'   for use with other functions of the \pkg{ThurstonianIRT} package.
#'
#' @examples
#' # load the data
#' data("triplets")
#'
#' # define the blocks of items
#' blocks <-
#'   set_block(c("i1", "i2", "i3"), traits = c("t1", "t2", "t3"),
#'           signs = c(1, 1, 1)) +
#'   set_block(c("i4", "i5", "i6"), traits = c("t1", "t2", "t3"),
#'             signs = c(-1, 1, 1)) +
#'   set_block(c("i7", "i8", "i9"), traits = c("t1", "t2", "t3"),
#'             signs = c(1, 1, -1)) +
#'   set_block(c("i10", "i11", "i12"), traits = c("t1", "t2", "t3"),
#'             signs = c(1, -1, 1))
#'
#' # generate the data to be understood by 'thurstonianIRT'
#' triplets_long <- make_TIRT_data(
#'   data = triplets, blocks = blocks, direction = "larger",
#'   format = "pairwise", family = "bernoulli", range = c(0, 1)
#' )
#'
#' \donttest{
#' # fit the data using Stan
#' fit <- fit_TIRT_stan(triplets_long, chains = 1)
#' print(fit)
#' predict(fit)
#' }
#'
#' @export
make_TIRT_data <- function(data, blocks, direction = c("larger", "smaller"),
                           format = c("ranks", "pairwise"),
                           family = "bernoulli", partial = FALSE,
                           range = c(0, 1)) {
  if (is.TIRTdata(data)) {
    return(data)
  }
  direction <- match.arg(direction)
  format <- match.arg(format)
  family <- check_family(family)
  stopifnot(is.TIRTblocks(blocks))
  blocks <- blocks$blocks
  data <- as.data.frame(data)
  npersons <- nrow(data)
  nblocks <- length(blocks)
  nitems_per_block <- length(blocks[[1]]$items)
  ncat <- 2L
  if (format == "ranks") {
    if (family != "bernoulli") {
      stop("Format 'ranks' requires family 'bernoulli'.")
    }
    partial <- as_one_logical(partial)
    # range is only meaningful for pairwise comparisons
    range <- c(0, 1)
  } else if (format == "pairwise") {
    if (length(range) != 2L) {
      stop("Argument 'range' must be of length 2.")
    }
    if (family == "cumulative") {
      ncat <- as.integer(range[2] - range[1] + 1)
    }
    # partial is only meaningful for ranked comparisons
    partial <- FALSE
  }
  ncomparisons <- (nitems_per_block * (nitems_per_block - 1)) / 2
  items_all <- ulapply(blocks, "[[", "items")
  if (any(duplicated(items_all))) {
    stop("Item variables in different blocks needs to have different names. ",
         "Use the 'names' argument in 'set_block' to equate item ",
         "parameters across blocks.")
  }
  nitems <- length(items_all)
  if (nitems != nitems_per_block * nblocks) {
    stop("All blocks should contain the same number of items.")
  }
  traits_all <- unique(ulapply(blocks, "[[", "traits"))
  ntraits <- length(traits_all)

  out <- tibble::tibble(
    person = rep(1:npersons, ncomparisons * nblocks),
    block = rep(1:nblocks, each = npersons * ncomparisons),
    comparison = rep(rep(1:ncomparisons, each = npersons), nblocks)
  )
  for (i in seq_len(nblocks)) {
    items <- blocks[[i]]$items
    item1 <- rep_comp(items, 1, nitems_per_block)
    item2 <- rep_comp(items, 2, nitems_per_block)
    traits <- blocks[[i]]$traits
    trait1 <- rep_comp(traits, 1, nitems_per_block)
    trait2 <- rep_comp(traits, 2, nitems_per_block)
    signs <- blocks[[i]]$signs
    sign1 <- rep_comp(signs, 1, nitems_per_block)
    sign2 <- rep_comp(signs, 2, nitems_per_block)

    rows <- out$block == i
    comparison <- out[rows, ]$comparison
    out[rows, "itemC"] <- comparison + (i - 1) * ncomparisons
    out[rows, "trait1"] <- trait1[comparison]
    out[rows, "trait2"] <- trait2[comparison]
    out[rows, "item1"] <- item1[comparison]
    out[rows, "item2"] <- item2[comparison]
    out[rows, "sign1"] <- sign1[comparison]
    out[rows, "sign2"] <- sign2[comparison]

    if (format == "ranks") {
      response1 <- unname(do.call(c, data[, item1, drop = FALSE]))
      response2 <- unname(do.call(c, data[, item2, drop = FALSE]))
      out[rows, "response1"] <- response1
      out[rows, "response2"] <- response2
      if (direction == "smaller") {
        response <- as.numeric(case_when(
          response1 < response2 ~ 1,
          response1 > response2 ~ 0,
          TRUE ~ NaN
        ))
      } else if (direction == "larger") {
        response <- as.numeric(case_when(
          response1 > response2 ~ 1,
          response1 < response2 ~ 0,
          TRUE ~ NaN
        ))
      }
    } else if (format == "pairwise") {
      item_comb_names <- paste0(item1, item2)
      response <- unname(do.call(c, data[, item_comb_names, drop = FALSE]))
      if (any(response < range[1] | response > range[2])) {
        stop("Responses must be within [", range[1], ", ", range[2], "].")
      }
      if (direction == "smaller") {
        # invert item responses
        response <- range[2] - response + range[1]
      }
      # the smallest value should be 0
      response <- response - range[1]
      if (family == "bernoulli") {
        # ensure dichotomous response format
        response <- ifelse(response > 0, 1, 0)
      }
    }
    out[rows, "response"] <- response
  }
  out$item1 <- factor(out$item1, levels = items_all)
  out$item2 <- factor(out$item2, levels = items_all)
  out$trait1 <- factor(out$trait1, levels = traits_all)
  out$trait2 <- factor(out$trait2, levels = traits_all)

  # check for partial comparisons
  if (format == "ranks") {
    is_nan <- is.nan(out$response)
    any_nan <- any(is_nan)
    if (any_nan) {
      if (!partial) {
        stop("Please set 'partial = TRUE' when using partial comparisons.")
      }
      out <- filter(out, !is_nan)
    } else {
      partial <- FALSE
    }
  }

  # check for items being used multiple times in the test
  item_names <- ulapply(blocks, "[[", "names")
  dupl_item_nums <- which(duplicated(item_names))
  dupl_item_names <- unique(item_names[dupl_item_nums])
  dupl_items <- named_list(dupl_item_names)
  for (i in seq_along(dupl_items)) {
    dupl_items[[i]] <- which(item_names == dupl_item_names[i])
  }
  # add attributes to the returned object
  structure(out,
    npersons = npersons,
    ntraits = ntraits,
    nblocks = nblocks,
    nitems = nitems,
    nitems_per_block = nitems_per_block,
    signs = ulapply(blocks, "[[", "signs"),
    dupl_items = dupl_items,
    traits = unique(ulapply(blocks, "[[", "traits")),
    format = format,
    family = family,
    partial = partial,
    range = range,
    ncat = ncat,
    class = c("TIRTdata", class(out))
  )
}

#' Prepare data for Thurstonian IRT models fitted with
#' lavaan or Mplus
#'
#' @inheritParams fit_TIRT_stan
#'
#' @return A \code{data.frame} ready to be passed to \pkg{lavaan}
#' or \pkg{Mplus}.
#'
#' @examples
#' # simulate some data
#' sdata <- sim_TIRT_data(
#'   npersons = 100,
#'   ntraits = 3,
#'   nblocks_per_trait = 4,
#'   gamma = 0,
#'   lambda = c(runif(6, 0.5, 1), runif(6, -1, -0.5)),
#'   Phi = diag(3)
#' )
#'
#' # create data ready for use in SEM software
#' sem_data <- make_sem_data(sdata)
#' head(sem_data)
#'
#' @import dplyr
#' @importFrom magrittr '%>%'
#' @export
make_sem_data <- function(data) {
  if (!is.TIRTdata(data)) {
    stop("'data' should be of class 'TIRTdata'. See ?make_TIRT_data")
  }
  data <- convert_factors(data)
  att <- attributes(data)
  npersons <- att[["npersons"]]
  ntraits <- att[["ntraits"]]
  nblocks <- att[["nblocks"]]
  ncols <- ntraits * (ntraits - 1) / 2 * nblocks
  data %>%
    mutate(itemC = paste0("i", .data$item1, "i", .data$item2)) %>%
    mutate(itemC = factor(.data$itemC, levels = unique(.data$itemC))) %>%
    select("person", "itemC", "response") %>%
    tidyr::spread(key = "itemC", value = "response") %>%
    select(-.data$person)
}

convert_factors <- function(data) {
  # data and code generating functions require
  # items and traits to be numeric
  stopifnot(is.TIRTdata(data))
  for (v in c("item1", "item2", "trait1", "trait2")) {
    data[[v]] <- as.integer(data[[v]])
  }
  data
}

#' Prepare blocks of items
#'
#' Prepare blocks of items and incorporate information
#' about which item belongs to which trait. A block
#' of items is a set of two or more items presented and answered together
#' by fully ranking them or selecting the most and/or least favorit
#' in a forced choice format. A whole test usually contains
#' several blocks and items may reappear in different blocks.
#'
#' @param items Names of item comparisons to be combined
#' into one block. Should correspond to variables in the data.
#' @param traits Names of the traits to which each item belongs
#' @param names Optional names of the items in the output.
#' Can be used to equate parameters of items across blocks,
#' if the same item was used in different blocks.
#' @param signs Expected signs of the item loadings (1 or -1).
#'
#' @examples
#' set_block(
#'   items = c("i1", "i2", "i3"),
#'   traits = c("A", "B", "C")
#' ) +
#' set_block(
#'   items = c("i4", "i5", "i6"),
#'   traits = c("A", "B", "C")
#' )
#'
#' # Support items i1 and i4 were the same so that they have the same parameters
#' set_block(
#'   items = c("i1", "i2", "i3"),
#'   traits = c("A", "B", "C"),
#'   names = c("item1", "item2", "item3")
#' ) +
#' set_block(
#'   items = c("i4", "i5", "i6"),
#'   traits = c("A", "B", "C"),
#'   names = c("item1", "item5", "item6")
#' )
#'
#' @seealso \code{\link{set_blocks_from_df}}
#' @export
set_block <- function(items, traits, names = items, signs = 1) {
  stopifnot(length(items) == length(traits))
  items <- as.character(items)
  traits <- as.character(traits)
  names <- as.character(names)
  if (length(signs) == 1L) {
    signs <- rep(signs, length(items))
  }
  stopifnot(length(items) == length(signs))
  signs <- sign(signs)
  out <- list(blocks = list(nlist(items, traits, names, signs)))
  structure(out, class = "TIRTblocks")
}

#' @rdname set_block
#' @export
empty_block <- function() {
  structure(list(blocks = list()), class = "TIRTblocks")
}

#' @export
"+.TIRTblocks" <- function(e1, e2) {
  stopifnot(is.TIRTblocks(e2))
  e1$blocks <- c(e1$blocks, e2$blocks)
  e1
}

#' Prepare blocks of items from a data frame
#'
#' Prepare blocks of items and incorporate information
#' about which item belongs to which trait from a pre-existing dataframe.
#' This is a wrapper function for \code{\link{set_block}}, eliminating the need
#' to manually set each item, trait, name and sign (loading) info per block.
#'
#' A block of items is a set of two or more items presented and answered
#' together by fully ranking them or selecting the most and/or least favorite
#' in a forced choice format. A whole test usually contains
#' several blocks and items may reappear in different blocks.
#'
#' @param data A \code{data.frame} containing all the required columns
#' (see the arguments below) to specify the item blocks.
#' @param blocks Name of column vector denoting the block each item
#' corresponds to. Each block must have an equal number of items.
#' @param items Name of column vector denoting items to be combined into
#' one block. Should correspond to variables in the data.
#' @param traits Names of column vector denoting the traits to which each
#' item belongs.
#' @param names Optional column vector of item names in the output.
#' Can be used to equate parameters of items across blocks,
#' if the same item was used in different blocks.
#' @param signs Name of column vector with expected signs of the
#' item loadings (1 or -1).
#'
#' @examples
#' block_info <- data.frame(
#'   block = rep(1:4, each = 3),
#'   items = c("i1", "i2", "i3", "i4", "i5", "i6",
#'             "i7", "i8", "i9", "i10", "i11", "i12"),
#'   traits = rep(c("t1", "t2", "t3"), times = 4),
#'   signs = c(1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, 1)
#' )
#'
#' blocks <- set_blocks_from_df(
#'   data = block_info,
#'   blocks = "block",
#'   items = "items",
#'   traits = "traits",
#'   signs = "signs"
#' )
#'
#' @seealso \code{\link{set_block}}
#' @export
set_blocks_from_df <- function(data, blocks = "block", items = "item",
                               traits = "trait", names = items,
                               signs = "sign") {
  # input checks
  data <- as.data.frame(data)
  blocks <- as_one_character(blocks)
  items <- as_one_character(items)
  traits <- as_one_character(traits)
  names <- as_one_character(names)
  signs <- as_one_character(signs)
  # check each block has the same number of items
  block_ids <- as.character(data[[blocks]])
  nitems_per_block <- table(as.character(data[[blocks]]))
  if (length(unique(nitems_per_block)) > 1L) {
    stop("All blocks should contain the same number of items.")
  }
  if (any(nitems_per_block < 2)) {
    stop("Blocks should contain more than one item.")
  }
  # save unique block_ids and number of blocks
  block_ids <- unique(block_ids)
  nblocks <- length(block_ids)
  # fill list with the set_block call for each block
  block_list <- list()
  for(i in seq_along(block_ids)) {
    sel <- data[[blocks]] == block_ids[i]
    block_list[[i]] <- set_block(
      items = data[sel, items], traits = data[sel, traits],
      names = data[sel, names], signs = data[sel, signs]
    )
  }
  # concatenate blocks
  Reduce("+", block_list)
}

is.TIRTdata <- function(x) {
  inherits(x, "TIRTdata")
}

is.TIRTblocks <- function(x) {
  inherits(x, "TIRTblocks")
}

rep_comp <- function(x, comp, nitems_per_block) {
  # generate comparisons of first vs. second items, traits etc.
  # Args:
  #   x: vector of names of items, traits, etc.
  #   comp: first or second comparison?
  stopifnot(comp %in% 1:2)
  if (comp == 1) {
    out <- rep(x[1:(nitems_per_block - 1)], (nitems_per_block - 1):1)
  } else {
    out <- ulapply(2:nitems_per_block, function(y) x[y:nitems_per_block])
  }
  out
}

check_family <- function(family, software = NULL) {
  options <- family_options(software)
  match.arg(family, options)
}

family_options <- function(software = NULL) {
  if (is.null(software)) {
    # TODO: the 'beta' family is implemented in Stan but still
    # needs to be understood theoretically before exporting it
    all_ops <- c("bernoulli", "cumulative", "gaussian")
    return(all_ops)
  }
  software <- match.arg(software, c("stan", "lavaan", "mplus"))
  if (software == "stan") {
    out <- c("bernoulli", "cumulative", "gaussian")
  } else if (software == "lavaan") {
    out <- c("bernoulli", "gaussian")
  } else if (software == "mplus") {
    out <- c("bernoulli", "gaussian")
  }
  out
}

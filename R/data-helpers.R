#' @export
make_TIRT_data <- function(data, blocks) {
  if (is.TIRTdata(data)) {
    return(data)
  }
  stopifnot(is.TIRTblocks(blocks))
  data <- as.data.frame(data)
  npersons <- nrow(data)
  nblocks <- length(blocks)
  nitems_per_block <- length(blocks[[1]]$items)
  ncomparisons <- (nitems_per_block * (nitems_per_block - 1)) / 2
  items_all <- ulapply(blocks, "[[", "items")
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
    item1 <- rep(
      items[1:(nitems_per_block - 1)], (nitems_per_block - 1):1
    )
    item2 <- unlist(lapply(
      2:nitems_per_block, function(x) items[x:nitems_per_block]
    ))
    traits <- blocks[[i]]$traits
    trait1 <- rep(
      traits[1:(nitems_per_block - 1)], (nitems_per_block - 1):1
    )
    trait2 <- unlist(lapply(
      2:nitems_per_block, function(x) traits[x:nitems_per_block]
    ))
    fblock <- (i - 1) * nitems_per_block
    comparison <- out[out$block == i, ]$comparison
    out[out$block == i, "itemC"] <- comparison + fblock
    out[out$block == i, "trait1"] <- trait1[comparison]
    out[out$block == i, "trait2"] <- trait2[comparison]
    out[out$block == i, "item1"] <- item1[comparison]
    out[out$block == i, "item2"] <- item2[comparison]
    resp_item1 <- unname(do.call(c, data[, item1]))
    resp_item2 <- unname(do.call(c, data[, item2]))
    out[out$block == i, "response"] <- as.numeric(resp_item1 < resp_item2)
  }
  out$item1 <- factor(out$item1, levels = items_all)
  out$item2 <- factor(out$item2, levels = items_all)
  out$trait1 <- factor(out$trait1, levels = traits_all)
  out$trait2 <- factor(out$trait2, levels = traits_all)
  
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
    npersons = npersons, ntraits = ntraits, nblocks = nblocks,
    nitems = nitems, nitems_per_block = nitems_per_block,
    signs = ulapply(blocks, "[[", "signs"), dupl_items = dupl_items,
    class = c("TIRTdata", class(out))
  )
}

#' @import dplyr
#' @importFrom magrittr '%>%'
#' @export
make_sem_data <- function(data, blocks = NULL) {
  if (!is.TIRTdata(data)) {
    data <- make_TIRT_data(data, blocks)
  }
  data <- convert_factors(data)
  att <- attributes(data)
  npersons <- att[["npersons"]]
  ntraits <- att[["ntraits"]] 
  nblocks <- att[["nblocks"]]
  ncols <- ntraits * (ntraits - 1) / 2 * nblocks
  data %>% 
    mutate(itemC = paste0("i", item1, "i", item2)) %>%
    mutate(itemC = factor(itemC, levels = unique(itemC))) %>%
    select(person, itemC, response) %>%
    tidyr::spread(key = itemC, value = response) %>%
    select(-person)
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
  structure(
    list(nlist(items, traits, names, signs)), 
    class = "TIRTblocks"
  )
}

#' @export 
empty_block <- function() {
  structure(list(), class = "TIRTblocks")
}

#' @export
"+.TIRTblocks" <- function(e1, e2) {
  stopifnot(is.TIRTblocks(e2))
  structure(c(e1, e2), class = "TIRTblocks")
}

is.TIRTdata <- function(x) {
  inherits(x, "TIRTdata")
}

is.TIRTblocks <- function(x) {
  inherits(x, "TIRTblocks")
}

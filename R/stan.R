#' @useDynLib thurstonianIRT, .registration = TRUE 
#' @import Rcpp
#' @import methods
#' @export
make_stan_data <- function(data, blocks = NULL) {
  if (!is.TIRTdata(data)) {
    data <- make_TIRT_data(data, blocks)
  }
  data <- convert_factors(data)
  att <- attributes(data)
  stan_data = list(
    Y = data$response,
    N = nrow(data),
    N_item = length(unique(c(data$item1, data$item2))),
    N_itemC = length(unique(data$itemC)),
    N_person = length(unique(data$person)),
    N_trait = length(unique(c(data$trait1, data$trait2))),
    J_item1 = as.numeric(data$item1),
    J_item2 = as.numeric(data$item2),
    J_itemC = as.numeric(data$itemC),
    J_person = as.numeric(data$person),
    J_trait1 = as.numeric(data$trait1),
    J_trait2 = as.numeric(data$trait2),
    J_item_pos = which(att[["signs"]] >= 0),
    J_item_neg = which(att[["signs"]] < 0)
  )
  nitems_per_block <- att[["nitems_per_block"]]
  nitems <- att[["nitems"]]
  stan_data$J_item_fix <- seq(1, nitems, nitems_per_block)
  stan_data$J_item_est <- setdiff(1:nitems, stan_data$J_item_fix)
  if (length(att$dupl_items)) {
    # force item parameters of the same item to be equal
    # this happens if the same items is applied in multiple blocks
    J_item_equal <- J_item_orig <- vector("list", length(att$dupl_items))
    for (i in seq_along(att$dupl_items)) {
      first <- att$dupl_items[[i]][1]
      dup <- att$dupl_items[[i]][-1]
      J_item_equal[[i]] <- dup
      J_item_orig[[i]] <- rep(first, length(dup))
    }
    stan_data$J_item_equal <- unlist(J_item_equal)
    stan_data$J_item_orig <- unlist(J_item_orig)
    # duplicated items should not be part of the other index variables
    stan_data$J_item_pos <- with(stan_data, setdiff(J_item_pos, J_item_equal))
    stan_data$J_item_neg <- with(stan_data, setdiff(J_item_neg, J_item_equal))
    stan_data$J_item_fix <- with(stan_data, setdiff(J_item_fix, J_item_equal))
    stan_data$J_item_est <- with(stan_data, setdiff(J_item_est, J_item_equal))
  } else { 
    stan_data$J_item_equal <-  stan_data$J_item_orig <- integer(0)
  }
  stan_data$N_item_pos = length(stan_data$J_item_pos)
  stan_data$N_item_neg = length(stan_data$J_item_neg)
  stan_data$N_item_fix = length(stan_data$J_item_fix)
  stan_data$N_item_est = length(stan_data$J_item_est)
  stan_data$N_item_equal <- length(stan_data$J_item_equal)
  stan_data$N_item_orig <- length(stan_data$J_item_orig)
  stan_data
}

#' @export
fit_TIRT_stan <- function(data, blocks = NULL, init = 0, ...) {
  data <- make_TIRT_data(data, blocks)
  stan_data <- make_stan_data(data)
  stan_pars = c("Cor_trait", "lambda", "psi", "gamma", "r", "eta")
  # TODO: fix problems when cores > 1
  rstan::sampling(
    stanmodels$thurstonian_irt_model, 
    data = stan_data, pars = stan_pars, 
    init = init, ...
  )
}

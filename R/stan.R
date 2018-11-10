#' Prepare data for Thurstonian IRT models fitted with Stan
#'
#' @inheritParams make_TIRT_data
#'
#' @return A list of data ready to be passed to \pkg{Stan}.
#'
#' @useDynLib thurstonianIRT, .registration = TRUE
#' @import Rcpp
#' @import methods
#' @export
make_stan_data <- function(data) {
  if (!is.TIRTdata(data)) {
    stop("'data' should be of class 'TIRTdata'. See ?make_TIRT_data")
  }
  data <- convert_factors(data)
  att <- attributes(data)
  out = list(
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
    J_item_neg = which(att[["signs"]] < 0),
    ncat = NCOL(data$gamma) + 1
  )

  # prepare family and response values
  family <- att$family
  options <- c("bernoulli", "cumulative", "beta", "normal")
  out$family <- as.numeric(factor(family, options))
  if (family %in% c("bernoulli", "cumulative")) {
    out$Yint <- data$response
    out$Yreal <- numeric(0)
  } else if (family %in% c("beta", "normal")) {
    out$Yint <- integer(0)
    out$Yreal <- data$response
  }

  nitems_per_block <- att[["nitems_per_block"]]
  nitems <- att[["nitems"]]
  out$J_item_fix <- seq(1, nitems, nitems_per_block)
  out$J_item_est <- setdiff(1:nitems, out$J_item_fix)
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
    out$J_item_equal <- unlist(J_item_equal)
    out$J_item_orig <- unlist(J_item_orig)
    # duplicated items should not be part of the other index variables
    out$J_item_pos <- with(out, setdiff(J_item_pos, J_item_equal))
    out$J_item_neg <- with(out, setdiff(J_item_neg, J_item_equal))
    out$J_item_fix <- with(out, setdiff(J_item_fix, J_item_equal))
    out$J_item_est <- with(out, setdiff(J_item_est, J_item_equal))
  } else {
    out$J_item_equal <-  out$J_item_orig <- integer(0)
  }
  out$N_item_pos = length(out$J_item_pos)
  out$N_item_neg = length(out$J_item_neg)
  out$N_item_fix = length(out$J_item_fix)
  out$N_item_est = length(out$J_item_est)
  out$N_item_equal <- length(out$J_item_equal)
  out$N_item_orig <- length(out$J_item_orig)
  out
}

#' Fit Thurstonian IRT models in Stan
#'
#' @param data An object of class \code{'TIRTdata'}. see
#' \code{\link{make_TIRT_data}} for documentation on how to create one.
#' @param init Initial values of the parameters.
#' Defaults to \code{0} as it proved to be most stable.
#' @param ... Further arguments passed to
#' \code{\link[rstan:sampling]{rstan::sampling}}.
#'
#' @return A \code{'TIRTfit'} object.
#'
#' @export
fit_TIRT_stan <- function(data, init = 0, ...) {
  stan_data <- make_stan_data(data)
  stan_pars = c(
    "Cor_trait", "lambda", "psi", "gamma",
    "gamma_ord", "disp", "r", "eta"
  )
  fit <- rstan::sampling(
    stanmodels$thurstonian_irt_model,
    data = stan_data, pars = stan_pars,
    init = init, ...
  )
  structure(nlist(fit, data), class = "TIRTfit")
}


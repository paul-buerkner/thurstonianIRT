#' Prepare data for Thurstonian IRT models fitted with Stan
#'
#' @inheritParams make_TIRT_data
#'
#' @return A list of data ready to be passed to \pkg{Stan}.
#'
#' #' @examples
#' # simulate some data
#' sim_data <- sim_TIRT_data(
#'   npersons = 100,
#'   ntraits = 3,
#'   nblocks_per_trait = 4,
#'   gamma = 0,
#'   lambda = c(runif(6, 0.5, 1), runif(6, -1, -0.5)),
#'   Phi = diag(3)
#' )
#'
#' # create data ready for use in Stan
#' stan_data <- make_stan_data(sim_data)
#' str(stan_data)
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
  # should be 2 for all non-ordinal families
  ncat <- if (!is.null(att$ncat)) att$ncat else 2L
  out = list(
    N = nrow(data),
    N_item = length(unique(c(data$item1, data$item2))),
    N_itemC = length(unique(data$itemC)),
    N_person = length(unique(data$person)),
    N_trait = length(unique(c(data$trait1, data$trait2))),
    J_item1 = as.array(as.numeric(data$item1)),
    J_item2 = as.array(as.numeric(data$item2)),
    J_itemC = as.array(as.numeric(data$itemC)),
    J_person = as.array(as.numeric(data$person)),
    J_trait1 = as.array(as.numeric(data$trait1)),
    J_trait2 = as.array(as.numeric(data$trait2)),
    J_item_pos = as.array(which(att[["signs"]] >= 0)),
    J_item_neg = as.array(which(att[["signs"]] < 0)),
    ncat = ncat
  )

  # prepare family and response values
  family <- check_family(att$family, "stan")
  options <- family_options("stan")
  out$family <- as.numeric(factor(family, options))
  if (family %in% c("bernoulli", "cumulative")) {
    out$Yint <- as.array(data$response)
    out$Yreal <- numeric(0)
  } else if (family %in% c("gaussian", "beta")) {
    out$Yint <- integer(0)
    out$Yreal <- as.array(data$response)
  }

  nitems_per_block <- att[["nitems_per_block"]]
  nitems <- att[["nitems"]]
  if (family %in% c("bernoulli", "cumulative", "beta")) {
    # fix first item uniqueness per block for identification
    # TODO: figure out if 'beta' really needs this
    out$J_item_fix <- as.array(seq(1, nitems, nitems_per_block))
  } else {
    out$J_item_fix <- integer(0)
  }
  out$J_item_est <- as.array(setdiff(1:nitems, out$J_item_fix))
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
    out$J_item_equal <- as.array(unlist(J_item_equal))
    out$J_item_orig <- as.array(unlist(J_item_orig))
    # duplicated items should not be part of the other index variables
    out$J_item_pos <- as.array(with(out, setdiff(J_item_pos, J_item_equal)))
    out$J_item_neg <- as.array(with(out, setdiff(J_item_neg, J_item_equal)))
    out$J_item_fix <- as.array(with(out, setdiff(J_item_fix, J_item_equal)))
    out$J_item_est <- as.array(with(out, setdiff(J_item_est, J_item_equal)))
  } else {
    out$J_item_equal <- out$J_item_orig <- integer(0)
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
fit_TIRT_stan <- function(data, init = 0, ...) {
  stan_data <- make_stan_data(data)
  stan_pars <- c(
    "Cor_trait", "lambda", "psi", "gamma",
    "gamma_ord", "disp", "r", "eta"
  )
  fit <- rstan::sampling(
    stanmodels$thurstonian_irt_model,
    data = stan_data, pars = stan_pars,
    init = init, ...
  )
  TIRTfit(fit, data)
}

# predict trait scores using Stan
predict_stan <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    out <- predict_stan_old_data(object, ...)
  } else {
    out <- predict_stan_new_data(object, newdata, ...)
  }
  out
}

# predict trait scores of existing persons using Stan
predict_stan_old_data <- function(object, ...) {
  fit <- object$fit
  traits <- attributes(object$data)$traits
  out <- as.data.frame(summary(fit, "eta")$summary)
  if (NROW(out)) {
    out <- out %>%
      tibble::rownames_to_column(var = "par") %>%
      rename(
        estimate = "mean", se = "sd",
        lower_ci = "2.5%", upper_ci = "97.5%"
      ) %>%
      select("par", "estimate", "se", "lower_ci", "upper_ci") %>%
      tidyr::extract(
        col = "par", into = c("par", "id", "trait"),
        regex = "(eta)\\[([[:digit:]]+),([[:digit:]]+)\\]"
      ) %>%
      mutate(
        id = as.integer(.data$id),
        trait = as.integer(.data$trait)
      ) %>%
      select(-.data$par) %>%
      mutate(trait = as.character(factor(.data$trait, labels = traits))) %>%
      arrange(.data$id)
  }
  as_tibble(out)
}

# predict trait scores of new persons using Stan
predict_stan_new_data <- function(object, newdata, inits = 0, ...) {
  # TODO: check 'newdata' for validity
  stan_data <- make_stan_data(newdata)

  # extract (medians of) item parameters and person hyperparameters
  # to be used as data in the predictions for new data
  par_dims <- object$fit@par_dims
  pars <- c("Cor_trait", "lambda", "psi", "r", "gamma", "gamma_ord", "disp")
  stan_data_pars <- named_list(pars)
  for (par in pars) {
    if (prod(par_dims[[par]]) > 0) {
      samples <- as.matrix(object$fit, par)
      stan_data_pars[[par]] <- apply(samples, 2, stats::median)
    } else {
      stan_data_pars[[par]] <- numeric(0)
    }
    dim(stan_data_pars[[par]]) <- par_dims[[par]]
  }
  # TODO: is there a more consistent approach to extract
  # the median cholesky factor of the trait correlation matrix?
  stan_data_pars$L_trait <- t(chol(stan_data_pars$Cor_trait))
  stan_data <- c(stan_data, stan_data_pars)

  # fit the model to obtain predictions for new persons
  fit <- rstan::sampling(
    stanmodels$thurstonian_irt_model_newdata,
    data = stan_data, pars = "eta",
    init = inits, ...
  )
  predict_stan_old_data(TIRTfit(fit, newdata), ...)
}

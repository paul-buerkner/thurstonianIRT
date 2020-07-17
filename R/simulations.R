#' Simulate Thurstonian IRT data
#'
#' @param npersons Number of persons.
#' @param ntraits Number of traits.
#' @param lambda Item factor loadings.
#' @param gamma Baseline attractiveness parameters of the
#'   first item versus the second item in the pairwise comparisons.
#'   Can be thought of as intercept parameters.
#' @param psi Optional item uniquenesses. If not provided,
#'   they will be computed as \code{psi = 1 - lambda^2} in which
#'   case lambda are taken to be the standardized factor loadings.
#' @param Phi Optional trait correlation matrix from which to sample
#'   person factor scores. Only used if \code{eta} is not provided.
#' @param eta Optional person factor scores. If provided, argument
#'   \code{Phi} will be ignored.
#' @param family Name of assumed the response distribution. Either
#'   \code{"bernoulli"}, \code{"cumulative"}, or \code{"gaussian"}.
#' @param nblocks_per_trait Number of blocks per trait.
#' @param nitems_per_block Number of items per block.
#' @param comb_blocks Indicates how to combine traits to blocks.
#'   \code{"fixed"} implies a simple non-random design that may combine
#'   certain traits which each other disproportionally often. We thus
#'   recommend to use a \code{"random"} block design (the default) that
#'   combines all traits with all other traits equally often on average.
#'
#' @return A \code{data.frame} of the same structure
#' as returned by \code{\link{make_TIRT_data}}. Parameter values
#' from which the data were simulated are stored as attributes
#' of the returned object.
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
#' # take a look at the data
#' head(sdata)
#' str(attributes(sdata))
#'
#' \donttest{
#' # fit a Thurstonian IRT model using lavaan
#' fit <- fit_TIRT_lavaan(sdata)
#' print(fit)
#' }
#'
#' @importFrom stats sd setNames rnorm rbeta
#' @importFrom rlang .data
#' @export
sim_TIRT_data <- function(npersons, ntraits, lambda, gamma,
                          psi = NULL, Phi = NULL, eta = NULL,
                          family = "bernoulli",
                          nblocks_per_trait = 5, nitems_per_block = 3,
                          comb_blocks = c("random", "fixed")) {
  # prepare data in long format to which responses may be added
  if ((ntraits * nblocks_per_trait) %% nitems_per_block != 0L) {
    stop("The number of items per block must divide ",
         "the number of total items.")
  }
  family <- check_family(family)
  comb_blocks <- match.arg(comb_blocks)
  nblocks <- ntraits * nblocks_per_trait / nitems_per_block
  nitems <- nitems_per_block * nblocks
  ncomparisons <- (nitems_per_block * (nitems_per_block - 1)) / 2
  data <- tibble::tibble(
    person = rep(1:npersons, ncomparisons * nblocks),
    block = rep(1:nblocks, each = npersons * ncomparisons),
    comparison = rep(rep(1:ncomparisons, each = npersons), nblocks)
  )
  # select traits for each block
  trait_combs <- make_trait_combs(
    ntraits, nblocks_per_trait, nitems_per_block,
    comb_blocks = comb_blocks
  )
  items_per_trait <- vector("list", ntraits)
  for (i in seq_len(nblocks)) {
    traits <- trait_combs[i, ]
    trait1 <- rep_comp(traits, 1, nitems_per_block)
    trait2 <- rep_comp(traits, 2, nitems_per_block)
    fblock <- (i - 1) * ncomparisons
    item1 <- match(trait1, traits) + fblock
    item2 <- match(trait2, traits) + fblock
    sign1 <- sign(lambda[item1])
    sign2 <- sign(lambda[item2])
    comparison <- data[data$block == i, ]$comparison
    data[data$block == i, "itemC"] <- comparison + fblock
    data[data$block == i, "trait1"] <- trait1[comparison]
    data[data$block == i, "trait2"] <- trait2[comparison]
    data[data$block == i, "item1"] <- item1[comparison]
    data[data$block == i, "item2"] <- item2[comparison]
    data[data$block == i, "sign1"] <- sign1[comparison]
    data[data$block == i, "sign2"] <- sign2[comparison]
    # save item numbers per trait
    for (t in unique(trait1)) {
      items_per_trait[[t]] <- union(
        items_per_trait[[t]], item1[match(t, trait1)]
      )
    }
    for (t in unique(trait2)) {
      items_per_trait[[t]] <- union(
        items_per_trait[[t]], item2[match(t, trait2)]
      )
    }
  }

  # prepare parameters
  if (is.null(eta)) {
    eta <- sim_eta(npersons, Phi)
  }
  if (length(gamma) == 1L) {
    gamma <- rep(gamma, ncomparisons * nblocks)
  }
  if (!is.list(lambda) && length(lambda) == 1L) {
    lambda <- rep(lambda, nitems)
  }
  if (is.null(psi)) {
    message("Computing standardized psi^2 as 1 - lambda^2")
    psi <- lambda2psi(lambda)
  } else if (!is.list(psi) && length(psi) == 1L) {
    psi <- rep(psi, nitems)
  }
  if (NROW(gamma) != ncomparisons * nblocks) {
    stop("gamma should contain ", ncomparisons * nblocks, " rows.")
  }
  if (sum(lengths(lambda)) != nitems) {
    stop("lambda should contain ", nitems, " values.")
  }
  if (sum(lengths(psi)) != nitems) {
    stop("psi should contain ", nitems, " values.")
  }
  dim_eta_exp <- c(length(unique(data$person)), ntraits)
  if (!is_equal(dim(eta), dim_eta_exp)) {
    stop("eta should be of dimension (", dim_eta_exp[1],
         ", ", dim_eta_exp[2], ").")
  }
  if (family == "cumulative") {
    stopifnot(NCOL(gamma) > 1L)
    data$gamma <- gamma[data$itemC, , drop = FALSE]
  } else {
    stopifnot(NCOL(gamma) == 1L)
    data$gamma <- as.vector(gamma)[data$itemC]
  }
  if (is.list(lambda)) {
    if (length(lambda) != ntraits) {
      stop("lambda should contain ", ntraits, " list elements.")
    }
    lambda_order <- order(unlist(items_per_trait))
    lambda <- unlist(lambda)[lambda_order]
  }
  data$lambda1 <- lambda[data$item1]
  data$lambda2 <- lambda[data$item2]
  if (is.list(psi)) {
    if (length(psi) != ntraits) {
      stop("psi should contain ", ntraits, " list elements.")
    }
    psi_order <- order(unlist(items_per_trait))
    psi <- unlist(psi)[psi_order]
  }
  data$psi1 <- psi[data$item1]
  data$psi2 <- psi[data$item2]
  for (p in seq_len(npersons)) {
    take <- data$person == p
    pdat <- data[take, ]
    data[take, "eta1"] <- eta[p, pdat$trait1]
    data[take, "eta2"] <- eta[p, pdat$trait2]
    # sample errors to make them item-person specific
    # but constant for each item per person
    errors <- rnorm(nitems, 0, psi)
    data[take, "error1"] <- errors[pdat$item1]
    data[take, "error2"] <- errors[pdat$item2]
  }

  data <- add_response(data, family = family)
  structure(data,
    npersons = npersons, ntraits = ntraits, nblocks = nblocks,
    nitems = nitems, nblocks_per_trait = nblocks_per_trait,
    nitems_per_block = nitems_per_block,
    signs = sign(lambda), lambda = lambda, psi = psi, eta = eta,
    traits = paste0("trait", seq_len(ntraits)),
    family = family, ncat = NCOL(data$gamma) + 1,
    class = c("TIRTdata", class(data))
  )
}

sim_eta <- function(npersons, Phi) {
  mu <- rep(0, nrow(Phi))
  mvtnorm::rmvnorm(npersons, mu, Phi)
}

add_response <- function(data, family) {
  # add columns related to the response
  # compute the latent mean 'mu'
  if (family %in% c("bernoulli", "gaussian", "beta")) {
    data$mu <- with(data, -gamma + lambda1 * eta1 - lambda2 * eta2)
  } else if (family %in% "cumulative") {
    # do not include 'gamma' here as it serves as the thresholds
    data$mu <- with(data, lambda1 * eta1 - lambda2 * eta2)
  }

  # deterministically include error sampled earlier
  mu_error <- with(data, mu + error1 - error2)
  data$error1 <- data$error2 <- NULL
  sum_psi <- with(data, sqrt(psi1^2 + psi2^2))

  # compute the actual response values
  if (family == "bernoulli") {
    data$response <- as.integer(mu_error >= 0)
  } else if (family == "cumulative") {
    data$response <- rep(NA, nrow(data))
    thres <- cbind(-Inf, data$gamma, Inf)
    for (i in seq_len(nrow(data))) {
      # 'cut' is not vectorized over 'breaks'
      data$response[i] <- cut(mu_error[i], breaks = thres[i, ])
    }
    # start counting at 0
    data$response <- as.integer(data$response) - 1
  } else if (family == "gaussian") {
    # do not use 'error' but sample directly from the latent distribution
    data$response <- rnorm(data$mu, sum_psi)
  } else if (family == "beta") {
    # mean parameterization of the beta distribution
    # do not use 'error' but sample directly from the latent distribution
    pr <- stats::pnorm(data$mu / sum_psi)
    data$response <- rbeta(length(pr), pr * data$disp, (1 - pr) * data$disp)
    # truncate distribution at the extremes
    data$response[data$response < 0.001] <- 0.001
    data$response[data$response > 0.999] <- 0.999
  }
  data
}

make_trait_combs <- function(ntraits, nblocks_per_trait, nitems_per_block,
                             comb_blocks = c("fixed", "random"),
                             maxtrys_outer = 20, maxtrys_inner = 1e6) {
  comb_blocks <- match.arg(comb_blocks)
  stopifnot((ntraits * nblocks_per_trait) %% nitems_per_block == 0L)
  if (comb_blocks == "fixed") {
    # use comb_blocks == "random" for a better balanced design
    traits <- rep(seq_len(ntraits), nblocks_per_trait)
    out <- matrix(traits, ncol = nitems_per_block, byrow = TRUE)
  } else if (comb_blocks == "random") {
    nblocks <- (ntraits * nblocks_per_trait) %/% nitems_per_block
    traits <- seq_len(ntraits)
    out <- replicate(nitems_per_block, traits, simplify = FALSE)
    out <- as.matrix(expand.grid(out))
    rownames(out) <- NULL
    remove <- rep(FALSE, nrow(out))
    for (i in seq_len(nrow(out))) {
      if (length(unique(out[i, ])) < ncol(out)) {
        remove[i] <- TRUE
      }
    }
    out <- out[!remove, ]
    possible_rows <- seq_len(nrow(out))
    nbpt_chosen <- rep(0, ntraits)

    .choose <- function(nblocks, maxtrys) {
      # finds suitable blocks
      chosen <- rep(NA, nblocks)
      i <- ntrys <- 1
      while (i <= nblocks && ntrys <= maxtrys) {
        ntrys <- ntrys + 1
        chosen[i] <- possible_rows[sample(seq_along(possible_rows), 1)]
        traits_chosen <- out[chosen[i], ]
        nbpt_chosen[traits_chosen] <- nbpt_chosen[traits_chosen] + 1
        valid <- max(nbpt_chosen) <= min(nbpt_chosen) + 1 &&
          !any(nbpt_chosen[traits_chosen] > nblocks_per_trait)
        if (valid) {
          possible_rows <- possible_rows[-chosen[i]]
          i <- i + 1
        } else {
          # revert number of blocks per trait chosen
          nbpt_chosen[traits_chosen] <- nbpt_chosen[traits_chosen] - 1
        }
      }
      return(chosen)
    }

    i <- 1
    chosen <- rep(NA, nblocks)
    while (anyNA(chosen) && i <= maxtrys_outer) {
      i <- i + 1
      chosen <- .choose(nblocks, maxtrys = maxtrys_inner)
    }
    if (anyNA(chosen)) {
      stop("Could not find a set of suitable blocks.")
    }
    out <- out[chosen, ]
  }
  out
}

lambda2psi <- function(lambda) {
  # according to Brown et al. 2011 for std lambda: psi^2 = 1 - lambda^2
  # psi is the SD and psi^2 is the variance
  # Ideas for lambda if unstandardized:
  # multiply std lambda with sqrt(2)
  # create simulated data based on std factor scores
  .lambda2psi <- function(x) {
    x <- as.numeric(x)
    if (any(abs(x) > 1)) {
      stop("standardized lambdas are expect to be between -1 and 1.")
    }
    sqrt(1 - x^2)
  }
  if (is.list(lambda)) {
    psi <- lapply(lambda, .lambda2psi)
  } else {
    psi <- .lambda2psi(lambda)
  }
  psi
}

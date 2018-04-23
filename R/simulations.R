#' Simulate Thurstonian IRT data
#' 
#' @param npersons Number of persons
#' @param ntraits Number of traits
#' @param gamma intercept parameters
#' @param lambda item factor loadings
#' @param psi optional item uniquenesses
#' @param Phi optional trait correlation matrix
#' @param eta optional person factor scores
#' @param nblocks_per_trait Number of blocks per trait
#' @param nitems_per_block Number of items per block
#' 
#' @return A \code{data.frame} of the same structure
#' as returned by \code{\link{make_TIRT_data}}.
#' 
#' @examples 
#' sdata <- sim_thurstonian_data(
#'   npersons = 100,
#'   ntraits = 3,
#'   nblocks_per_trait = 4,
#'   gamma = 0,
#'   lambda = c(runif(6, 0.5, 1), runif(6, -1, -0.5)),
#'   Phi = diag(3)
#' )
#' head(sdata)
#' 
#' @importFrom stats sd setNames
#' @importFrom rlang .data
#' @export
sim_thurstonian_data <- function(npersons, ntraits, gamma, lambda, 
                                 psi = NULL, Phi = NULL, eta = NULL, 
                                 nblocks_per_trait = 5, nitems_per_block = 3) {
  # prepare data in long format to which responses may be added
  if ((ntraits * nblocks_per_trait) %% nitems_per_block != 0L) {
    stop("The number of items per block must divide ", 
         "the number of total items.")
  }
  nblocks <- ntraits * nblocks_per_trait / nitems_per_block
  nitems <- nitems_per_block * nblocks
  ncomparisons <- (nitems_per_block * (nitems_per_block - 1)) / 2
  data <- tibble::tibble(
    person = rep(1:npersons, ncomparisons * nblocks),
    block = rep(1:nblocks, each = npersons * ncomparisons),
    comparison = rep(rep(1:ncomparisons, each = npersons), nblocks)
  )
  # select traits for each block
  trait_combs <- make_trait_combs(ntraits, nblocks_per_trait, nitems_per_block)
  items_per_trait <- vector("list", ntraits)
  for (i in seq_len(nblocks)) {
    traits <- trait_combs[i, ]
    trait1 <- rep(
      traits[1:(nitems_per_block - 1)], (nitems_per_block - 1):1
    )
    trait2 <- unlist(lapply(
      2:nitems_per_block, function(x) traits[x:nitems_per_block]
    ))
    fblock <- (i - 1) * nitems_per_block
    item1 <- match(trait1, traits) + fblock
    item2 <- match(trait2, traits) + fblock
    comparison <- data[data$block == i, ]$comparison 
    data[data$block == i, "itemC"] <- comparison + fblock
    data[data$block == i, "trait1"] <- trait1[comparison] 
    data[data$block == i, "trait2"] <- trait2[comparison] 
    data[data$block == i, "item1"] <- item1[comparison] 
    data[data$block == i, "item2"] <- item2[comparison]
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
  if (length(gamma) != ncomparisons * nblocks) {
    stop("gamma should contain ", ncomparisons * nblocks, " values.")
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
  data$gamma <- gamma[data$comparison]
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
  }
  
  data$prob <- prob_response(data)
  data$response <- sim_response(data)
  structure(data, 
    npersons = npersons, ntraits = ntraits, nblocks = nblocks,
    nitems = nitems, nblocks_per_trait = nblocks_per_trait, 
    nitems_per_block = nitems_per_block, 
    signs = sign(lambda), lambda = lambda, psi = psi, eta = eta, 
    traits = paste0("trait", seq_len(ntraits)),
    class = c("TIRTdata", class(data))
  )
}

sim_eta <- function(npersons, Phi) {
  mu <- rep(0, nrow(Phi))
  mnormt::rmnorm(npersons, mu, Phi)
}

sim_response <- function(data) {
  prob <- prob_response(data)
  stats::rbinom(length(prob), size = 1, prob = prob)
}

prob_response <- function(data) {
  z <- with(data, 
    (- gamma + lambda1 * eta1 - lambda2 * eta2) / sqrt(psi1^2 + psi2^2)
  )
  stats::pnorm(z)
}

make_trait_combs <- function(ntraits, nblocks_per_trait, nitems_per_block) {
  # TODO: improve design balance
  stopifnot((ntraits * nblocks_per_trait) %% nitems_per_block == 0L)
  traits <- rep(seq_len(ntraits), nblocks_per_trait)
  matrix(traits, ncol = nitems_per_block, byrow = TRUE)
}

lambda2psi <- function(lambda) {
  # according to Brown et al. 2011 for std lambda
  # Ideas for lambda if unstandardized:
  # multiply std lambda with sqrt(2)
  # create simulated data based on std factor scores
  .lambda2psi <- function(x) {
    x <- as.numeric(x)
    if (any(abs(x) > 1)) {
      stop("standardized lambdas are expect to be between -1 and 1.")
    }
    1 - x^2
  }
  if (is.list(lambda)) {
    psi <- lapply(lambda, .lambda2psi)
  } else {
    psi <- .lambda2psi(lambda)
  }
  psi
}

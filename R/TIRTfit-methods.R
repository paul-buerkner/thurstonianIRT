TIRTfit <- function(fit, data) {
  version <- utils::packageVersion("thurstonianIRT")
  structure(nlist(fit, data, version), class = "TIRTfit")
}

is.TIRTfit <- function(x) {
  inherits(x, "TIRTfit")
}

#' @export
print.TIRTfit <- function(x, ...) {
  print(x$fit, ...)
}

#' @method summary TIRTfit
#' @importMethodsFrom lavaan summary
#' @export
summary.TIRTfit <- function(object, ...) {
  summary(object$fit, ...)
}

#' @export
predict.TIRTfit <- function(object, ...) {
  fit <- object$fit
  traits <- attributes(object$data)$traits
  if (inherits(fit, "stanfit")) {
    # post process Stan objects
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
  } else {
    # post process Mplus objects
    if (inherits(fit, "mplusObjectTIRT")) {
      out <- fit$results[["trait_scores"]]
      if (is.null(out)) {
        # for backwards compatibility with version < 0.9.3
        out <- fit$results[["savedata"]]
      }
      out <- as.data.frame(out)
      if (NROW(out)) {
        ntraits <- ncol(out)
        out <- out %>%
          tidyr::gather("trait", "estimate", everything()) %>%
          mutate(id = rep(seq_len(n() / ntraits), ntraits)) %>%
          arrange(.data$id) %>%
          select("id", "trait", "estimate")
      }
      se <- as.data.frame(fit$results[["trait_scores_se"]])
      if (NROW(se)) {
        ntraits <- ncol(se)
        se <- se %>%
          tidyr::gather("trait", "se", everything()) %>%
          mutate(id = rep(seq_len(n() / ntraits), ntraits)) %>%
          arrange(.data$id) %>%
          select("id", "trait", "se")
      }
      out <- out %>% inner_join(se, by = c("id", "trait"))
    } else if (inherits(fit, "lavaan")) {
      # post process lavaan objects
      out <- as.data.frame(lavaan::lavPredict(fit, ...))
      if (NROW(out)) {
        ntraits <- ncol(out)
        out <- out %>%
          tidyr::gather("trait", "estimate", everything()) %>%
          mutate(id = rep(seq_len(n() / ntraits), ntraits)) %>%
          arrange(.data$id) %>%
          select("id", "trait", "estimate")
      }
    }
  }
  as_tibble(out)
}

#' Extract corrected goodness of fit statistics
#'
#' By default \pkg{lavaan} will return a value for degrees of
#' freedom that ignores redundancies amongst the estimated model
#' thresholds. This function corrects the degrees of freedom, and
#' then recalculates the associated chi-square test statistic
#' p-value and root mean square error of approximation (RMSEA).
#'
#' Note this function is currently only implemented for \pkg{lavaan}.
#'
#' @param object A \code{TIRTfit} object.
#' @param ... currently unused.
#'
#' @return A vector containing the chi-square value, adjusted degrees of
#'   freedom, p-value, and RMSEA.
#'
#' @examples
#' # load the data
#' data("triplets")
#'
#' # define the blocks of items
#' blocks <-
#'   set_block(c("i1", "i2", "i3"), traits = c("t1", "t2", "t3"),
#'             signs = c(1, 1, 1)) +
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
#' # fit the data using lavaan
#' fit <- fit_TIRT_lavaan(triplets_long)
#' gof(fit)
#'
#' @export
gof.TIRTfit <- function(object, ...) {
  if (!inherits(object$fit, "lavaan")) {
    stop("gof.TIRTfit currently only works for lavaan fitted TIRT models.")
  }

  # Extract chi_sq, N, and unadjusted DF
  chi_sq <- object$fit@test$scaled.shifted$stat
  N <- length(unique(object$data$person))
  df <- object$fit@test$scaled.shifted$df

  # Get number of items per block to calculate redundancies
  blocks <- unique(object$data$block)
  redundancies <- rep(NA, length(blocks))
  for (i in seq_along(blocks)) {
    n_items <- nrow(unique(subset(object$data, block == blocks[i], select = itemC)))
    redundancies[i] <- n_items * (n_items - 1) * (n_items - 2) / 6
  }

  # Adjust the DF, p-value, and recalculate the RMSEA
  df <- df - sum(redundancies)
  p_val <- 1 - pchisq(chi_sq, df)
  RMSEA <- ifelse(df > chi_sq, 0, sqrt((chi_sq - df)/(df * (N - 1))))
  gof <- c(Chi_Sq = chi_sq, df = df, p_val = p_val, RMSEA = RMSEA)
  gof
}

#' @rdname gof.TIRTfit
#' @export
gof <- function(object, ...) {
  UseMethod("gof")
}

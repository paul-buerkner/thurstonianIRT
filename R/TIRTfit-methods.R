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

#' Predict trait scores of Thurstonian IRT models
#'
#' @param object An object of class \code{TIRTfit}.
#' @param newdata Optional \code{TIRTdata} object (created via
#'   \code{\link{make_TIRT_data}}) containing data of new persons
#'   for which trait scores should be predicted based on the fitted
#'   model. If \code{NULL} (the default), trait scores are predicted
#'   for the persons whose data was used to originally fit the model.
#' @param ... Further arguments passed to the underlying methods.
#'
#' @details When predicting trait scores of new persons (via \code{newdata}),
#'   posterior medians of item parameters are used for predictions. This implies
#'   that the uncertainty in the new trait scores is underestimated as the
#'   uncertainty in the (posterior distribution of) item parameters is ignored.
#'
#' @return A data frame with predicted trait scores.
#'
#' @export
predict.TIRTfit <- function(object, newdata = NULL, ...) {
  if (inherits(object$fit, "stanfit")) {
    out <- predict_stan(object, newdata = newdata, ...)
  } else if (inherits(object$fit, "mplusObjectTIRT")) {
    out <- predict_mplus(object, newdata = newdata, ...)
  } else if (inherits(object$fit, "lavaan")) {
    out <- predict_lavaan(object, newdata = newdata, ...)
  }
  out
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

  # Extract N, chi_sq, and unadjusted DF
  N <- length(unique(object$data$person))
  chi_sq <- object$fit@test$scaled.shifted$stat
  if (is.null(chi_sq)) {
    chi_sq <- NA
  }
  df <- object$fit@test$scaled.shifted$df
  if (is.null(df)) {
    df <- NA
  }

  # Get number of items per block to calculate redundancies
  blocks <- unique(object$data$block)
  redundancies <- rep(NA, length(blocks))
  for (i in seq_along(blocks)) {
    sel_items <- unique(object$data$itemC[object$data$block == blocks[i]])
    n_items <- length(sel_items)
    redundancies[i] <- n_items * (n_items - 1) * (n_items - 2) / 6
  }

  # Adjust the DF, p-value, and recalculate the RMSEA
  df <- df - sum(redundancies)
  p_val <- 1 - stats::pchisq(chi_sq, df)
  RMSEA <- ifelse(df > chi_sq, 0, sqrt((chi_sq - df)/(df * (N - 1))))
  gof <- c(chi_sq = chi_sq, df = df, p_val = p_val, RMSEA = RMSEA)
  gof
}

#' @rdname gof.TIRTfit
#' @export
gof <- function(object, ...) {
  UseMethod("gof")
}

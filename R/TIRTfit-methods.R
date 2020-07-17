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

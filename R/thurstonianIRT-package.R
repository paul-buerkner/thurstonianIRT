#' The 'thurstonianIRT' package.
#'
#' @description This package fits Thurstonian Item Response Theory (IRT) models
#'   using 'Stan', 'lavaan', or 'Mplus'. To bring your data into the right
#'   format, use the \code{\link{make_TIRT_data}} function. Models can then be
#'   fitted via \code{\link{fit_TIRT_stan}}, \code{\link{fit_TIRT_lavaan}}, or
#'   \code{\link{fit_TIRT_mplus}} depending on the desired model fitting engine.
#'   Data from Thurstonian IRT models can be simulated via
#'   \code{\link{sim_TIRT_data}}.
#'
#' @docType package
#' @name thurstonianIRT-package
#' @aliases thurstonianIRT
#' @useDynLib thurstonianIRT, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#'
#' @references
#' Brown, A., & Maydeu-Olivares, A. (2011). Item response modeling of
#' forced-choice questionnaires. Educational and Psychological Measurement,
#' 71(3), 460-502. doi:10.1177/0013164410375112
#'
#' Bürkner P. C., Schulte N., & Holling H. (2019). On the Statistical and
#' Practical Limitations of Thurstonian IRT Models. Educational and
#' Psychological Measurement. doi:10.1177/0013164419832063
#'
NULL

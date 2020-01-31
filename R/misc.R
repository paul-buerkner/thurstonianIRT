bias <- function(x, y) {
  colMeans(x - y)
}

RMSE <- function(x, y) {
  sqrt(colMeans((x - y)^2))
}

#' Extract corrected goodness of fit statistics
#'
#' By default \code{lavaan} will return a value for degrees of
#' freedom that ignores redundancies amongst the estimated model
#' thresholds. This function corrects the degrees of freedom, and
#' then recalculates the associated chi-square test statistic
#' p-value and root mean square error of approximation (RMSEA).
#'
#' Note this function is currently only implemented for \code{lavaan}.
#'
#' @param TIRTfit a TIRTfit object
#'
#' @return A vector containing the chi-square value, adjusted degrees of freedom, p-value, and RMSEA.
#'
#' @examples
#' gof_TIRT(fit)
#'
#' @export
gof_TIRT <- function(TIRTfit){
  # Stop conditions: Must be TIRTfit obj & fitted in lavaan
  if(inherits(TIRTfit, "TIRTfit") == FALSE){
    stop("The first argument must be a lavaan-fitted TIRT model.")}
  if(inherits(TIRTfit$fit, "lavaan") == FALSE){
    stop("gof_TIRT currently only works for lavaan fitted TIRT models.")}

  # Extract chi_sq, N, and unadjusted DF
  chi_sq <- TIRTfit$fit@test$scaled.shifted$stat
  N <- length(unique(TIRTfit$data$person))
  df <- TIRTfit$fit@test$scaled.shifted$df

  # Get number of items per block to calculate redundancies
  block_no <- unique(TIRTfit$data$block)
  redundancies <- vector()
  for(y in seq(block_no)){
    n_items <- nrow(unique(subset(TIRTfit$data, block == y, select = itemC)))
    redundancies[y] <- n_items * (n_items - 1) * (n_items - 2) / 6
  }

  # Adjust the DF, p-value, and recalculate the RMSEA
  df <- df - sum(redundancies)
  p_val <- 1 - pchisq(chi_sq, df)
  ifelse(df > chi_sq,
         RMSEA <- 0,
         RMSEA <- sqrt((chi_sq - df)/(df * (N - 1))) )
  gof <- c(Chi_Sq = chi_sq, df = df, p_val = p_val, RMSEA = RMSEA)
  return(gof)
}

scale2 <- function(x, center = FALSE) {
  if (center) {
    means <- colMeans(x)
    x <- sweep(x, 2, means, "-")
  }
  sds <- apply(x, 2, sd)
  sweep(x, 2, sds, "/")
}

collapse <- function(..., sep = "") {
  # wrapper for paste with collapse = ""
  paste(..., sep = sep, collapse = "")
}

is_wholenumber <- function(x, tol = .Machine$double.eps) {
  # check if x is a whole number (integer)
  if (!is.numeric(x)) {
    out <- FALSE
  } else {
    out <- abs(x - round(x)) < tol
  }
  out
}

is_equal <- function(x, y, ...) {
  isTRUE(all.equal(x, y, ...))
}

#' Set up Correlation Matrices
#'
#' @param cors vector of unique correlations
#' @param dim Dimension of the correlation matrix
#' @param dimnames Optional dimnames of the correlation matrix
#'
#' @return A correlation \code{matrix} of dimension \code{dim}.
#'
#' @examples
#' cor_matrix(c(0.2, 0.3, 0.5), dim = 3)
#'
#' @export
cor_matrix <- function(cors, dim, dimnames = NULL) {
  out <- diag(dim)
  out[lower.tri(out)] <- cors
  out[upper.tri(out)] <- t(out)[upper.tri(out)]
  if (!is.null(dimnames)) {
    dimnames(out) <- list(dimnames, dimnames)
  }
  out
}

collapse_lines <- function(...) {
  dots <- c(...)
  paste0(dots, collapse = "\n")
}

nlist <- function(...) {
  # create a named list using object names
  m <- match.call()
  dots <- list(...)
  no_names <- is.null(names(dots))
  has_name <- if (no_names) FALSE
  else nzchar(names(dots))
  if (all(has_name)) return(dots)
  nms <- as.character(m)[-1]
  if (no_names) {
    names(dots) <- nms
  } else {
    names(dots)[!has_name] <- nms[!has_name]
  }
  dots
}

ulapply <- function(X, FUN, ...) {
  # short for unlist(lapply(.))
  unlist(lapply(X = X, FUN = FUN, ...))
}

named_list <- function(names, values = NULL) {
  # initialize a named list
  # Args:
  #   names: names of the elements
  #   values: values of the elements
  if (!is.null(values)) {
    if (length(values) <= 1L) {
      values <- replicate(length(names), values)
    }
    values <- as.list(values)
    stopifnot(length(values) == length(names))
  } else {
    values <- vector("list", length(names))
  }
  setNames(values, names)
}

# coerce 'x' to a single character string
as_one_character <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- as.character(x)
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- deparse_combine(s, max_char = 100L)
    stop2("Cannot coerce ", s, " to a single character value.")
  }
  x
}

# coerce 'x' to a signle number value
as_one_numeric <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- substr(deparse_combine(s), 1L, 100L)
    stop("Cannot coerce ", s, " to a single numeric value.")
  }
  x
}

# coerce 'x' to TRUE or FALSE if possible
as_one_logical <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- as.logical(x)
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- substr(deparse_combine(s), 1L, 100L)
    stop("Cannot coerce ", s, " to a single logical value.")
  }
  x
}

# combine deparse lines into one string
deparse_combine <- function(x, max_char = 100) {
  out <- collapse(deparse(x))
  if (isTRUE(max_char > 0)) {
    out <- substr(out, 1, max_char)
  }
  out
}

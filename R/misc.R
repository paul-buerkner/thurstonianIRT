bias <- function(x, y) {
  colMeans(x - y)
}

RMSE <- function(x, y) {
  sqrt(colMeans((x - y)^2))
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


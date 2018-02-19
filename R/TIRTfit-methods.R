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
    out <- summary(fit, "eta")$summary %>%
      as.data.frame() %>%
      mutate(par = rownames(.data)) %>%
      rename(
        estimate = "mean", se = "sd",
        lower_ci = "2.5%", upper_ci = "97.5%"
      ) %>%
      select("par", "estimate", "se", "lower_ci", "upper_ci") %>%
      tidyr::extract(
        col = "par", into = c("par", "id", "trait"),
        regex = "(eta)\\[([[:digit:]]+),([[:digit:]]+)\\]"
      ) %>%
      select(-.data$par) %>%
      mutate(trait = as.character(factor(.data$trait, labels = traits))) %>%
      arrange(.data$id) 
  } else {
    if (inherits(fit, "mplusObjectTIRT")) {
      out <- fit$results[["savedata"]]
    } else if (inherits(fit, "lavaan")) {
      out <- lavaan::lavPredict(fit, ...)
    }
    ntraits <- ncol(out)
    out <- as.data.frame(out) %>%
      tidyr::gather("trait", "estimate", names(.data)) %>%
      mutate(id = rep(seq_len(n() / ntraits), ntraits)) %>%
      arrange(.data$id) %>%
      select("id", "trait", "estimate")
  }
  as_tibble(out)
}

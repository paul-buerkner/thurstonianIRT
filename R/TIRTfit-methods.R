is.TIRTfit <- function(x) {
  inherits(x, "TIRTfit")
}

#' @export
print.TIRTfit <- function(x, ...) {
  print(x$fit, ...)
}

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
      mutate(par = rownames(.)) %>%
      rename(
        estimate = mean, se = sd,
        lower_ci = "2.5%", upper_ci = "97.5%"
      ) %>%
      select(par, estimate, se, lower_ci, upper_ci) %>%
      tidyr::extract(
        col = "par", into = c("par", "id", "trait"),
        regex = "(eta)\\[([[:digit:]]+),([[:digit:]]+)\\]"
      ) %>%
      select(-par) %>%
      mutate(trait = as.character(factor(trait, labels = traits))) %>%
      arrange(id) 
  } else {
    if (inherits(fit, "mplusObjectTIRT")) {
      out <- fit$results[["savedata"]]
    } else if (inherits(fit, "lavaan")) {
      out <- lavaan::lavPredict(fit, ...)
    }
    ntraits <- ncol(out)
    out <- as.data.frame(out) %>%
      tidyr::gather("trait", "estimate", names(.)) %>%
      mutate(id = rep(seq_len(n() / ntraits), ntraits)) %>%
      arrange(id) %>%
      select(id, trait, estimate)
  }
  as_tibble(out)
}

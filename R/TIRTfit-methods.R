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


# Get tidy df of model parameters from lavaam
TIRT_parameters <- function(TIRTfit){

  # Get information about the blocks
  block_data <- TIRTfit$data %>%  # Reduce df to one row per item
    distinct(item1, item2, .keep_all = TRUE)
  nitems <- max(as.numeric(block_data$item2))  # Item 2 always has last item
  nblocks <- max(as.numeric(block_data$block))
  nitems_per_block <- nitems / nblocks

  # Get TIRT parameters from lavaan
  parameters <- lavaan::parameterEstimates(TIRTfit$fit)

  # Get the signs to adjust loadings
  # Need to stack both item and sign columns and get unique rows
  signs <- tibble(item = with(block_data, c(item1, item2)),
                  sign = with(block_data, c(sign1, sign2))) %>%
    dplyr::distinct() %>%
    dplyr::arrange(item) %>%
    dplyr::pull(sign)

  # Extract loadings and adjust by the sign
  loadings <- parameters %>%
    select(label, lhs, est, se) %>%
    slice(seq(1, nitems * (nitems_per_block - 1),
              by = nitems_per_block - 1)) %>%
    mutate(label = as.numeric(gsub("\\D", "", label))) %>%
    arrange(label) %>%
    mutate(est = abs(est) * signs) %>%
    rename(trait = lhs, item_number = label)

  # Get pairwise variance parameters (psi_i^2 + psi_k^2)
  # This is the diagonal of the error covariance matrix
  vars <- parameters %>%
    select(label, est, se) %>%
    slice(grep(".P.", label)) %>%  #" "^[^P]*P[^P]*$" for uniquenesses
    mutate(item1 = as.numeric(sub("P", "", sub("P[0-9]+$", "", label))),
           item2 = as.numeric(sub(".*P", "", label))) %>%
    arrange(item1, item2) %>%
    select(item1, item2, est, se)

  # Get item thresholds
  thresholds <- parameters %>%
    filter(rhs == "t1") %>%
    select(lhs, est, se) %>%
    mutate(item1 = as.numeric(sub("i", "", sub("i[0-9]+$", "", lhs))),
           item2 = as.numeric(sub(".*i", "", lhs))) %>%
    arrange(item1, item2) %>%
    select(item1, item2, est, se)

  return(
    list(
      loadings = loadings,
      variances = vars,
      thresholds = thresholds
    )
  )
}

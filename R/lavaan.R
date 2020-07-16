#' Generate lavaan code for Thurstonian IRT models
#'
#' @inheritParams fit_TIRT_stan
#'
#' @return A character string of lavaan code
#' for a Thurstonian IRT model.
#'
#' @examples
#' lambdas <- c(runif(6, 0.5, 1), runif(6, -1, -0.5))
#' sim_data <- sim_TIRT_data(
#'   npersons = 100,
#'   ntraits = 3,
#'   nblocks_per_trait = 4,
#'   gamma = 0,
#'   lambda = lambdas,
#'   Phi = diag(3)
#' )
#' cat(make_lavaan_code(sim_data))
#'
#' @export
make_lavaan_code <- function(data) {
  if (!is.TIRTdata(data)) {
    stop("'data' should be of class 'TIRTdata'. See ?make_TIRT_data")
  }
  data <- convert_factors(data)
  data <- filter(data, .data$person == unique(.data$person)[1])
  att <- attributes(data)
  family <- check_family(att$family, "lavaan")
  nitems <- att[["nitems"]]
  nitems_per_block <- att[["nitems_per_block"]]
  ntraits <- att[["ntraits"]]
  traits <- seq_len(ntraits)
  if (isTRUE(att[["partial"]])) {
    stop("Cannot yet handle partial comparisons when using lavaan.")
  }

  # define factor loadings (lambda)
  lav_loadings <- vector("list", ntraits)
  for (i in traits) {
    for (n in seq_len(nrow(data))) {
      if (data$trait1[n] == i) {
        lav_loadings[[i]] <- c(lav_loadings[[i]], with(data,
           paste0("start(1) * L", item1[n], " * i", item1[n], "i", item2[n])
        ))
      } else if (data$trait2[n] == i) {
        lav_loadings[[i]] <- c(lav_loadings[[i]], with(data,
           paste0("start(-1) * L", item2[n], "n * i", item1[n], "i", item2[n])
        ))
      }
    }
    lav_loadings[[i]] <- paste0(
      "trait", i, " =~ ",
      paste0(lav_loadings[[i]], collapse = " + ")
    )
  }
  lav_loadings <- collapse(unlist(lav_loadings), "\n")

  # fix factor varianaces to 1
  lav_fix_factor_variances <-
    collapse("trait", traits, " ~~ 1 * trait", traits, "\n")

  # factor correlations
  lav_factor_correlations <- collapse(
    sapply(1:(ntraits - 1),
      function(i) paste0(
        "trait", i, " ~~ ",
        paste0("trait", (i + 1):ntraits, collapse = " + "),
        "\n"
     )
    )
  )

  # fix factor loadings of the same item to the same value
  lav_fix_factor_loadings <- ""
  items_both_dir <- which(1:nitems %in% data$item1 & 1:nitems %in% data$item2)
  if (length(items_both_dir)) {
    lav_fix_factor_loadings <- collapse(
      "L", items_both_dir, " == -L", items_both_dir, "n\n"
    )
  }

  # declare uniquenesses (psi)
  lav_uniqueness <- with(data, collapse(
    "i", item1, "i", item2,
    " ~~ P", item1, "P", item2,
    " * i", item1, "i", item2, "\n"
  ))

  # correlated uniqunesses
  lav_cor_uniqueness <- ""
  for (n in 1:(nrow(data) - 1)) {
    for (m in (n + 1):nrow(data)) {
      pos_psi1 <- with(data, item1[n] == item1[m])
      pos_psi2 <- with(data, item2[n] == item2[m])
      neg_psi <- with(data, item2[n] == item1[m])
      if (pos_psi1) {
        lav_cor_uniqueness <- with(data,
          paste0(lav_cor_uniqueness,
            "i", item1[n], "i", item2[n], " ~~ ",
            "start(1) * P", item1[n],
            " * i", item1[m], "i", item2[m], "\n"
          )
        )
      } else if (pos_psi2) {
        lav_cor_uniqueness <- with(data,
          paste0(lav_cor_uniqueness,
            "i", item1[n], "i", item2[n], " ~~ ",
            "start(1) * P", item2[n],
            " * i", item1[m], "i", item2[m], "\n"
          )
        )
      } else if (neg_psi) {
        lav_cor_uniqueness <- with(data,
          paste0(lav_cor_uniqueness,
            "i", item1[n], "i", item2[n], " ~~ ",
            "start(-1) * P", item2[n], "n",
            " * i", item1[m], "i", item2[m], "\n"
          )
        )
      }
    }
  }

  # pair's uniqueness is equal to sum of 2 utility uniqunesses
  lav_equal_uniqueness <- ""
  if (nitems_per_block > 2) {
    psi_item1 <- paste0("P", data$item1)
    psi_item2 <- paste0("P", data$item2)
    neg_psi1 <- sapply(paste0(" ", psi_item1, "n "), grepl, lav_cor_uniqueness)
    neg_psi2 <- sapply(paste0(" ", psi_item2, "n "), grepl, lav_cor_uniqueness)
    lav_equal_uniqueness <- with(data, collapse(
      psi_item1, psi_item2, " == ",
      ifelse(neg_psi1, paste0(" - ", psi_item1, "n"), psi_item1),
      ifelse(neg_psi2, paste0(" - ", psi_item2, "n"), paste0(" + ", psi_item2)),
      "\n"
    ))
  }

  # fix certain uniquenesses for identification
  lav_fix_uniqueness <- ""
  if (family %in% "bernoulli") {
    if (nitems_per_block > 2) {
      # fix one uniqueness per block for identification
      lav_fix_uniqueness <- collapse(
        "P", seq(1, nitems, nitems_per_block), " == 1\n"
      )
    } else {
      # fix all uniquenesses for identification
      psi_item1 <- paste0("P", data$item1)
      psi_item2 <- paste0("P", data$item2)
      lav_fix_uniqueness <- collapse(psi_item1, psi_item2, " == 1\n")
    }
  }

  # force item parameters of the same item to be equal
  # this happens if the same items is applied in multiple blocks
  lav_equal_items <- ""
  for (i in seq_along(att$dupl_items)) {
    first <- att$dupl_items[[i]][1]
    dup <- att$dupl_items[[i]][-1]
    lav_equal_items <- paste0(lav_equal_items,
      collapse("L", first, " == L", dup, "\n"),
      collapse("P", first, " == P", dup, "\n")
    )
  }

  # combine all lavaan code snippets
  collapse_lines(
    "# factor loadings (lambda)",
    lav_loadings,
    "# fix factor variances to 1",
    lav_fix_factor_variances,
    "# factor correlations",
    lav_factor_correlations,
    "# fix factor loadings of the same item to the same value",
    lav_fix_factor_loadings,
    "# declare uniquenesses (psi)",
    lav_uniqueness,
    "# correlated uniqunesses",
    lav_cor_uniqueness,
    "# pair's uniqueness is equal to sum of 2 utility uniqunesses",
    lav_equal_uniqueness,
    "# fix certain uniquenesses for identification",
    lav_fix_uniqueness,
    "# force item parameters of the same item to be equal",
    lav_equal_items
  )
}

#' Fit Thurstonian IRT models in lavaan
#'
#' @inheritParams fit_TIRT_stan
#' @param estimator Name of the estimator that should be used.
#'   See \code{\link[lavaan:lavOptions]{lavOptions}}.
#' @param ... Further arguments passed to
#'   \code{\link[lavaan:lavaan]{lavaan}}.
#'
#' @return A \code{'TIRTfit'} object.
#'
#' @examples
#' # load the data
#' data("triplets")
#'
#' # define the blocks of items
#' blocks <-
#'   set_block(c("i1", "i2", "i3"), traits = c("t1", "t2", "t3"),
#'           signs = c(1, 1, 1)) +
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
#' \donttest{
#' # fit the data using lavaan
#' fit <- fit_TIRT_lavaan(triplets_long)
#' print(fit)
#' predict(fit)
#' }
#'
#' @export
fit_TIRT_lavaan <- function(data, estimator = "ULSMV", ...) {
  lavaan_data <- make_sem_data(data)
  lavaan_model <- make_lavaan_code(data)

  att <- attributes(data)
  family <- check_family(att$family, "lavaan")
  if (family %in% "bernoulli") {
    ordered <- names(lavaan_data)
  } else if (family %in% "gaussian") {
    ordered <- NULL
  }
  fit <- lavaan::lavaan(
    lavaan_model, data = lavaan_data, ordered = ordered,
    auto.fix.first = FALSE, auto.th = TRUE,
    parameterization = "theta", estimator = estimator,
    ...
  )
  TIRTfit(fit, data)
}

# predict trait scores using lavaan
predict_lavaan <- function(object, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    # TODO: check 'newdata' for validity
    newdata <- make_sem_data(newdata)
  }
  fit <- object$fit
  traits <- attributes(object$data)$traits
  out <- as.data.frame(lavaan::lavPredict(fit, newdata = newdata, ...))
  if (NROW(out)) {
    ntraits <- ncol(out)
    out <- out %>%
      tidyr::gather("trait", "estimate", everything()) %>%
      mutate(id = rep(seq_len(n() / ntraits), ntraits)) %>%
      arrange(.data$id) %>%
      select("id", "trait", "estimate")
  }
  as_tibble(out)
}

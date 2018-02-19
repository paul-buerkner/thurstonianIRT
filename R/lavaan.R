#' Generate lavaan code for Thurstonian IRT models
#' 
#' @inheritParams make_TIRT_data
#' 
#' @return A character string of lavaan code
#' for a Thurstonian IRT model.
#' 
#' @examples 
#' lambdas <- c(runif(6, 0.5, 1), runif(6, -1, -0.5))
#' sdata <- sim_thurstonian_data(
#'   npersons = 100,
#'   ntraits = 3,
#'   nblocks_per_trait = 4,
#'   gamma = 0,
#'   lambda = lambdas,
#'   Phi = diag(3)
#' )
#' cat(make_lavaan_code(sdata))
#' 
#' @export
make_lavaan_code <- function(data, blocks = NULL) {
  if (!is.TIRTdata(data)) {
    data <- make_TIRT_data(data, blocks)
  }
  data <- convert_factors(data)
  data <- filter(data, .data$person == unique(.data$person)[1])
  att <- attributes(data)
  nitems <- att[["nitems"]] 
  nitems_per_block <- att[["nitems_per_block"]]
  ntraits <- att[["ntraits"]]
  traits <- seq_len(ntraits)
  
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
        paste0("trait", (i+1):ntraits, collapse = " + "),
        "\n"
     )
    )
  )
  
  # fix factor loadings of the same item to the same value
  items_both_dir <- which(
    1:nitems %in% data$item1 & 1:nitems %in% data$item2
  )
  lav_fix_factor_loadings <- collapse(
    "L", items_both_dir, " == -L", items_both_dir, "n\n"
  )
  
  # declare uniquenesses (psi)
  lav_uniqueness <- with(data, collapse(
    "i", item1, "i", item2, 
    " ~~ P", item1, "P", item2, 
    " * i", item1, "i", item2, "\n"
  ))
  
  # correlated uniqunesses
  lav_cor_uniqueness <- ""
  for (n in 1:(nrow(data) - 1)) {
    for (m in (n+1):nrow(data)) {
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
  
  # fix one uniqueness per block for identification
  lav_fix_uniqueness <- collapse(
    "P", seq(1, nitems, nitems_per_block), " == 1\n"
  )
  
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
    "# fix one uniqueness per block for identification",
    lav_fix_uniqueness,
    "# force item parameters of the same item to be equal",
    lav_equal_items
  )
}

#' Fit Thurstonian IRT models in lavaan
#' 
#' @inheritParams make_TIRT_data
#' @param estimator Name of the estimator that should be used. 
#'   See \code{\link[lavaan:lavOptions]{lavOptions}}.
#' @param ... Further arguments passed to 
#'   \code{\link[lavaan:lavaan]{lavaan}}.
#' 
#' @return A \code{TIRTfit} object.
#' 
#' @export
fit_TIRT_lavaan <- function(data, blocks = NULL, estimator = "ULSMV", ...) {
  data <- make_TIRT_data(data, blocks)
  lavaan_data <- make_sem_data(data)
  lavaan_model <- make_lavaan_code(data)
  fit <- lavaan::lavaan(
    lavaan_model, data = lavaan_data, ordered = names(lavaan_data),
    auto.fix.first = FALSE, auto.th = TRUE,
    parameterization = "theta", estimator = estimator, 
    ...
  )
  structure(nlist(fit, data), class = "TIRTfit")
}

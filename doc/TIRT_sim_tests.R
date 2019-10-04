params <-
list(EVAL = TRUE)

## ---- SETTINGS-knitr, include=FALSE-----------------------------------------------------
stopifnot(require(knitr))
options(width = 90)
opts_chunk$set(
  comment = NA,
  message = FALSE,
  warning = FALSE,
  eval = if (isTRUE(exists("params"))) params$EVAL else FALSE,
  dev = "png",
  dpi = 150,
  fig.asp = 0.8,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center"
)

## ---------------------------------------------------------------------------------------
library(thurstonianIRT)
library(dplyr)
library(tidyr)

## ---------------------------------------------------------------------------------------
set.seed(1234)
npersons <- 500
ntraits <- 5
nitems_per_block <- 3
nblocks_per_trait <- 9
nblocks <- ntraits * nblocks_per_trait / nitems_per_block
nitems <- ntraits * nblocks_per_trait
ncomparisons <- (nitems_per_block * (nitems_per_block - 1)) / 2 * nblocks

## ---------------------------------------------------------------------------------------
set.seed(1234)
lambda <- runif(nitems, 0.65, 0.96)
signs <- c(rep(1, ceiling(nitems / 2)), rep(-1, floor(nitems / 2)))
lambda <- lambda * signs[sample(seq_len(nitems))]
gamma <- runif(nitems, -1, 1)
Phi <- diag(5)

## ---------------------------------------------------------------------------------------
sdata <- sim_TIRT_data(
  npersons = npersons, 
  ntraits = ntraits, 
  nitems_per_block = nitems_per_block,
  nblocks_per_trait = nblocks_per_trait,
  gamma = gamma,
  lambda = lambda,
  Phi = Phi
)

## ---- results="hide"--------------------------------------------------------------------
fit_stan <- fit_TIRT_stan(sdata, chains = 1, iter = 1000, warmup = 500)
fit_lavaan <- fit_TIRT_lavaan(sdata)
fit_mplus <- fit_TIRT_mplus(sdata)

## ---------------------------------------------------------------------------------------
eta <- as_tibble(as.data.frame(attributes(sdata)$eta))
names(eta) <- paste0("trait", 1:ncol(eta))
true_scores <- eta %>%
  mutate(id = 1:n()) %>%
  gather(key = "trait", value = "truth", -id)
true_summaries <- true_scores %>%
  group_by(trait) %>%
  summarise(true_mean = mean(truth), true_sd = sd(truth))

pred <- predict(fit_stan) %>% 
  bind_rows(predict(fit_lavaan), predict(fit_mplus), .id = "source") %>%
  mutate(
    source = as.character(factor(
      source, levels = 1:3, labels = c("stan", "lavaan", "mplus")
    )),
    trait = tolower(trait)
  ) %>%
  inner_join(true_scores, by = c("id", "trait"))

pred <- pred %>%
  inner_join(
    pred %>%
      group_by(trait, source) %>%
      summarise(cor_est_truth = cor(estimate, truth)), 
    by = c("trait", "source")
  ) %>%
  mutate(
    sign = sign(cor_est_truth),
    estimate = ifelse(sign %in% -1, -estimate, estimate)
  ) %>%
  inner_join(true_summaries, by = "trait") %>%
  group_by(trait, source) %>%
  mutate(
    est_mean = mean(estimate),
    est_sd = sd(estimate)
  ) %>%
  ungroup() %>%
  mutate(
    ztruth = (truth - true_mean) / true_sd,
    zestimate = (estimate - est_mean) / est_sd
  )

## ---------------------------------------------------------------------------------------
res <- pred %>%
  group_by(trait, source) %>%
  summarise(rel = cor(estimate, truth)^2)

res

## ---- include = FALSE-------------------------------------------------------------------
testthat::expect_true(all(res$rel > 0.84))

## ---------------------------------------------------------------------------------------
cor_matrix <- pred %>%
  mutate(
    # ensure correct ordering of traits
    SC = paste0(source, "_", trait),
    SC = factor(SC, levels = unique(SC))
  ) %>%
  select(id, SC, estimate) %>%
  spread(key = "SC", value = "estimate") %>%
  bind_cols(eta, .) %>%
  select(-id) %>%
  cor()

## ---------------------------------------------------------------------------------------
trait1 <- paste0(c("stan", "lavaan", "mplus"), "_trait1")
round(cor_matrix[trait1, trait1], 2)

## ---- include = FALSE-------------------------------------------------------------------
for (i in 1:ntraits) {
  trait_cols <- paste0(c("stan", "lavaan", "mplus"), "_trait", i)
  testthat::expect_true(all(cor_matrix[trait_cols, trait_cols] > 0.97))
}


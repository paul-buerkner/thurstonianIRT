functions {
  /* cumulative log-PDF for a single response
   * assumes the latent variable to be normal
   * Args:
   *   y: response category
   *   mu: linear predictor
   *   thres: ordinal thresholds
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real cumulative_Phi_lpmf(int y, real mu, vector thres) {
     int ncat = num_elements(thres) + 1;
     real p;
     if (y == 0) {
       p = Phi(thres[1] - mu);
     } else if (y == ncat - 1) {
       p = 1 - Phi(thres[ncat - 1] - mu);
     } else {
       p = Phi(thres[y + 1] - mu) - Phi(thres[y] - mu);
     }
     return log(p);
   }
}
data {
  // model type
  int<lower=1,upper=4> family;
  // family == 1: bernoulli
  // family == 2: cumulative
  // family == 3: gaussian
  // family == 4: beta
  int<lower=1> N;  // total number of observations
  // response variable
  int Yint[family == 1 || family == 2 ? N : 0];
  real Yreal[family == 3 || family == 4 ? N : 0];
  int<lower=1> N_item;
  int<lower=1> N_itemC;  // item pairs
  int<lower=1> N_person;
  int<lower=1> N_trait;
  int<lower=0> N_item_fix;
  int<lower=0> N_item_est;
  // indices over N
  int<lower=1> J_item1[N];
  int<lower=1> J_item2[N];
  int<lower=1> J_itemC[N];
  int<lower=1> J_person[N];
  int<lower=1> J_trait1[N];
  int<lower=1> J_trait2[N];
  int<lower=1> J_item_fix[N_item_fix];
  int<lower=1> J_item_est[N_item_est];
  // indicate inverted items
  int<lower=0> N_item_pos;
  int<lower=0> N_item_neg;
  int<lower=1> J_item_pos[N_item_pos];
  int<lower=1> J_item_neg[N_item_neg];
  // indicate items used in multiple blocks
  int<lower=0> N_item_equal;
  int<lower=1> J_item_equal[N_item_equal];
  int<lower=1> J_item_orig[N_item_equal];
  // number of response categories in ordinal models
  // should be set to 2 for other families
  int<lower=2> ncat;
}
transformed data {
  vector<lower=0>[N_item_fix] psi_fix;  // fixed item SDs
  psi_fix = rep_vector(1.0, N_item_fix);
}
parameters {
  // item thresholds depend on the family
  vector[family == 1 || family == 3 || family == 4 ? N_itemC : 0] gamma;
  ordered[ncat - 1] gamma_ord[family == 2 ? N_itemC : 0];
  vector<lower=0>[N_item_pos] lambda_pos;  // item loadings
  vector<upper=0>[N_item_neg] lambda_neg;  // item loadings
  vector<lower=0>[N_item_est] psi_est;  // estimated item SDs
  // unscaled group level effects of traits within persons
  matrix[N_trait, N_person] z_trait;
  // cholesky factor of correlation matrix of traits
  cholesky_factor_corr[N_trait] L_trait;
  vector[N_item] z;  // unscaled random effects
  // dispersion parameter of the beta family
  // TODO: make comparison specific?
  vector<lower=0>[family == 4 ? 1 : 0] disp;
}
transformed parameters {
  // latent traits per person
  matrix[N_person, N_trait] eta;
  vector[N_item] lambda;  // item loadings
  vector<lower=0>[N_item] psi;  // item SDs
  vector[N_item] r;  // item random effects
  eta = (L_trait * z_trait)';
  psi[J_item_fix] = psi_fix;
  psi[J_item_est] = psi_est;
  psi[J_item_equal] = psi[J_item_orig];
  lambda[J_item_pos] = lambda_pos;
  lambda[J_item_neg] = lambda_neg;
  lambda[J_item_equal] = lambda[J_item_orig];
  // accounts for the correlations within the same block
  r = psi .* z;
}
model {
  vector[N] mu;
  vector[N] sum_psi;
  for (n in 1:N) {
    // compute linear predictor
    mu[n] = r[J_item1[n]] - r[J_item2[n]] +
      lambda[J_item1[n]] * eta[J_person[n], J_trait1[n]] -
      lambda[J_item2[n]] * eta[J_person[n], J_trait2[n]];
    // compute item-comparison standard deviations
    sum_psi[n] = sqrt(psi[J_item1[n]]^2 + psi[J_item2[n]]^2);
  }
  // likelihood contributions
  if (family == 1) {
    // bernoulli models
    for (n in 1:N) {
      // use - gamma for consistency with Brown et al. 2011
      mu[n] = Phi((mu[n] - gamma[J_itemC[n]]) / sum_psi[n]);
    }
    Yint ~ bernoulli(mu);
  } else if (family == 2) {
    // cumulative models
    for (n in 1:N) {
      // scale quantities
      vector[ncat - 1] thres = gamma_ord[J_itemC[n]] / sum_psi[n];
      mu[n] /= sum_psi[n];
      // likelihood contribution
      Yint[n] ~ cumulative_Phi(mu[n], thres);
    }
  } else if (family == 3) {
    // gaussian models
    for (n in 1:N) {
      mu[n] = (mu[n] - gamma[J_itemC[n]]) / sum_psi[n];
    }
    Yreal ~ normal(mu, 1);
  } else if (family == 4) {
    // beta models
    for (n in 1:N) {
      mu[n] = Phi((mu[n] - gamma[J_itemC[n]]) / sum_psi[n]);
    }
    Yreal ~ beta(mu * disp[1], (1 - mu) * disp[1]);
  }
  // prior specifications
  if (family == 2) {
    for (i in 1:N_itemC) {
      gamma_ord[i] ~ normal(0, 3);
    }
  } else {
    gamma ~ normal(0, 3);
  }
  lambda_pos ~ normal(1, 0.5);
  lambda_neg ~ normal(-1, 0.5);
  psi_est ~ normal(1, 0.3);
  L_trait ~ lkj_corr_cholesky(1);
  to_vector(z_trait) ~ normal(0, 1);
  z ~ normal(0, 1);
}
generated quantities {
  // Cor_trait is Phi
  corr_matrix[N_trait] Cor_trait;
  Cor_trait = multiply_lower_tri_self_transpose(L_trait);
}

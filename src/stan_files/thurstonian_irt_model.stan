functions {
#include /chunks/functions-cumulative_Phi_lpmf.stan
}
data {
#include /chunks/data-shared.stan
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
  // TODO: make comparison-specific?
  vector<lower=0>[family == 4 ? 1 : 0] disp;
}
transformed parameters {
  // latent traits per person
  matrix[N_person, N_trait] eta;
  vector[N_item] lambda;  // item loadings
  vector<lower=0>[N_item] psi;  // item SDs
  vector[N_item] r;  // item random effects (residuals)
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
#include /chunks/model-likelihood.stan
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

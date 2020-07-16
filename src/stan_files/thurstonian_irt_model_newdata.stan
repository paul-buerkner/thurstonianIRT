functions {
#include /chunks/functions-cumulative_Phi_lpmf.stan
}
data {
#include /chunks/data-shared.stan
  // fix item parameters and person hyperparameters
  // for predictions of parameters for new persons
  vector[N_item] lambda;  // item loadings
  vector<lower=0>[N_item] psi;  // item SDs
  vector[N_item] r;  // item random effects (residuals)
  // item thresholds depend on the family
  vector[family == 1 || family == 3 || family == 4 ? N_itemC : 0] gamma;
  ordered[ncat - 1] gamma_ord[family == 2 ? N_itemC : 0];
  // dispersion parameter of the beta family
  vector<lower=0>[family == 4 ? 1 : 0] disp;
  // cholesky factor of correlation matrix of traits
  cholesky_factor_corr[N_trait] L_trait;
}
transformed data {
}
parameters {
  // unscaled group level effects of traits within persons
  matrix[N_trait, N_person] z_trait;
}
transformed parameters {
  // latent traits per person
  matrix[N_person, N_trait] eta;
  eta = (L_trait * z_trait)';
}
model {
#include /chunks/model-likelihood.stan
  to_vector(z_trait) ~ normal(0, 1);
}
generated quantities {
}

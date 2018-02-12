functions { 
} 
data { 
  int<lower=1> N;  // total number of observations 
  int Y[N];  // response variable 
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
  // fix lambda
  // vector[N_item] lambda;
} 
transformed data {
  vector<lower=0>[N_item_fix] psi_fix;  // fixed item SDs
  psi_fix = rep_vector(1.0, N_item_fix);
} 
parameters { 
  vector[N_itemC] gamma;  // item thresholds
  vector<lower=0>[N_item_pos] lambda_pos;  // item loadings
  vector<upper=0>[N_item_neg] lambda_neg;  // item loadings
  vector<lower=0>[N_item_est] psi_est;  // estimated item SDs
  // unscaled group level effects of traits within persons
  matrix[N_trait, N_person] z_trait;
  // cholesky factor of correlation matrix of traits
  cholesky_factor_corr[N_trait] L_trait; 
  vector[N_item] z;  // unscaled random effects
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
  for (n in 1:N) { 
    // define linear predictor
    mu[n] = - gamma[J_itemC[n]] + 
      lambda[J_item1[n]] * eta[J_person[n], J_trait1[n]] -
      lambda[J_item2[n]] * eta[J_person[n], J_trait2[n]] + 
      r[J_item1[n]] - r[J_item2[n]];
    // scale with item variances
    mu[n] = mu[n] / sqrt(psi[J_item1[n]]^2 + psi[J_item2[n]]^2);
    // transform to probability space
    mu[n] = Phi(mu[n]);
  } 
  // prior specifications 
  gamma ~ normal(0, 3);
  lambda_pos ~ normal(1, 0.5);
  lambda_neg ~ normal(-1, 0.5);
  psi_est ~ normal(1, 0.3);
  L_trait ~ lkj_corr_cholesky(1); 
  to_vector(z_trait) ~ normal(0, 1);
  z ~ normal(0, 1);
  // likelihood contribution 
  Y ~ bernoulli(mu);
} 
generated quantities { 
  // Cor_trait is Phi
  corr_matrix[N_trait] Cor_trait;
  Cor_trait = multiply_lower_tri_self_transpose(L_trait); 
}

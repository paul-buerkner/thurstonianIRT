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
      mu[n] = mu[n] - gamma[J_itemC[n]];
    }
    // SD = 1 is arbitrary as the residuals r are already part of mu
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

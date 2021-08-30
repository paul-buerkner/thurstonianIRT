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

// BayGMST_R_v1.0: AR(1) structure in the temperature and proxy equations
//
// Ported unchanged from BayGMST_v1.0.stan (see inst/legacy-scripts/) as part
// of the BayGMST 0.1.0 package restructuring. Loaded at runtime via
// instantiate::stan_package_model(name = "baygmst", package = "BayGMST").

data {
  int<lower=1> NT;          // total time points
  int<lower=0> NT_obs;      // number of instrumental observations
  int<lower=0> NT_mis;      // number of missing years to reconstruct

  // indices to map observed and missing data to the full time vector
  array[NT_obs] int<lower=1, upper=NT> idx_obs;
  array[NT_mis] int<lower=1, upper=NT> idx_mis;

  // climate forcings
  vector[NT] G; // greenhouse gases
  vector[NT] S; // solar irradiance
  vector[NT] V; // volcanic activity

  // data
  vector[NT_obs] y_obs; // observed instrumental temperatures (T)
  vector[NT] z;         // proxy records (R)
}
parameters {
  real alpha0;
  real alpha1;
  real<lower=0> sigma_z; // noise in proxy equation (epsilon)
  real<lower=-1, upper=1> phi_R; // AR(1) coeff. for proxy

  real beta0;
  real betaG;
  real betaS;
  real betaV;

  real<lower=0> sigma_y; // noise in temperature equation (eta/nu)
  real<lower=-1, upper=1> phi_T; // AR(1) coeff. for temperature

  // declare the missing temperature values as parameters to be estimated
  vector[NT_mis] y_mis;
}
transformed parameters {
  // construct the full Temperature vector T (called 'y' here)
  // this merges the parameters (missing) and data (observed)
  vector[NT] y;
  y[idx_obs] = y_obs;
  y[idx_mis] = y_mis; // performing inference on the missing values of T

  // calculate the deterministic mean based on forcings
  vector[NT] mu_forcing;
  mu_forcing = beta0 + betaG * G + betaS * S + betaV * V;
}
model {
  // #### priors ####
  // #### priors ####
  alpha0 ~ normal(0, 0.5);
  alpha1 ~ normal(0, 0.5);
  sigma_z ~ exponential(1);
  phi_R ~ normal(0, 0.3);

  beta0 ~ normal(0, 0.5);
  betaG ~ normal(0, 0.5);
  betaS ~ normal(0, 0.5);
  betaV ~ normal(0, 0.5);
  sigma_y ~ exponential(1);
  phi_T ~ normal(0, 0.3);

  // a. temperature model (lagged dependent variable) ####
  // equation: y_t = phi_T * y_{t-1} + mu_forcing_t + error

  // t=1: we can't use t-1, so we estimate it near the forcing mean (should be good enough)
  y[1] ~ normal(mu_forcing[1], sigma_y);

  // t=2 to NT: The current T depends on the previous T
  y[2:NT] ~ normal(
    phi_T * y[1:NT-1] + mu_forcing[2:NT],
    sigma_y
  );

  // b. proxy model (lagged dependent variable) ####
  vector[NT] mu_proxy = alpha0 + alpha1 * y;

  // t=1
  z[1] ~ normal(mu_proxy[1], sigma_z);

  // t=2 to NT: the current proxy depends on the previous proxy
  z[2:NT] ~ normal(
    phi_R * z[1:NT-1] + mu_proxy[2:NT],
    sigma_z
  );
}
generated quantities {
  vector[NT] mu_y;              // fitted mean for y_t
  vector[NT_obs] y_ins_fitted;  // proxy-only estimate of the instrumental years (posterior predictive checks)

  // y mean (conditional on realized y[t-1])
  mu_y[1] = mu_forcing[1];
  for (t in 2:NT)
    mu_y[t] = phi_T * y[t-1] + mu_forcing[t];

  // y_ins_fitted: draw the instrumental-period temperatures from their
  // distribution given the proxies (not y_obs), the parameters, and the
  // reconstructed temperature just before the instrumental period. Kalman
  // filter forward, then backward sampling (forward-filtering
  // backward-sampling), so each year uses proxies from the whole period
  // rather than only earlier years. Contributed by Tyler Bagwell
  // (2026-11-24); it replaces a forward-only filter that inflated the
  // credible band at the start of the instrumental period.
  //
  // The filter runs over every year from the first to the last observed one
  // (L years), so gaps inside the instrumental record are handled; if the
  // record starts at t = 1 the filter starts from the prior on y[1].
  {
    int k0 = idx_obs[1];
    int L = idx_obs[NT_obs] - k0 + 1;
    vector[L] a;   // predicted mean
    vector[L] Pp;  // predicted variance
    vector[L] m;   // filtered mean
    vector[L] P;   // filtered variance
    vector[L] y_s; // backward-sampled path
    real m_prev = 0;
    real P_prev = 0;
    if (k0 > 1)
      m_prev = y[k0 - 1];

    for (j in 1:L) {
      int k = k0 + j - 1;
      real r;
      real K;
      if (k == 1) {
        a[j]  = mu_forcing[1];
        Pp[j] = square(sigma_y);
        r = z[1] - alpha0;
      } else {
        a[j]  = phi_T * m_prev + mu_forcing[k];
        Pp[j] = square(phi_T) * P_prev + square(sigma_y);
        r = z[k] - phi_R * z[k - 1] - alpha0;
      }
      K = Pp[j] * alpha1 / (square(alpha1) * Pp[j] + square(sigma_z));
      m[j] = a[j] + K * (r - alpha1 * a[j]);
      P[j] = (1 - K * alpha1) * Pp[j];
      m_prev = m[j];
      P_prev = P[j];
    }

    y_s[L] = normal_rng(m[L], sqrt(P[L]));
    for (s in 1:(L - 1)) {
      int j = L - s;
      real J = phi_T * P[j] / Pp[j + 1];
      // P[j] * (1 - J * phi_T), written in a form that stays positive
      y_s[j] = normal_rng(m[j] + J * (y_s[j + 1] - a[j + 1]),
                          sqrt(P[j] * square(sigma_y) / Pp[j + 1]));
    }

    for (t in 1:NT_obs)
      y_ins_fitted[t] = y_s[idx_obs[t] - k0 + 1];
  }
}

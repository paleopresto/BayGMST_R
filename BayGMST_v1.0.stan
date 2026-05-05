// BayGMST_v0.3: Same as BayGMST_v0.1 + AR(1) structure in the temperature and proxy equations

data {
  int<lower=1> NT;          // total time points
  int<lower=0> NT_obs;      // number of instrumental observations
  int<lower=0> NT_mis;      // number of missing years to reconstruct
  int<lower=0> NT_prj;      // number of years to project into the future
  
  // indices to map observed and missing data to the full time vector
  array[NT_obs] int<lower=1, upper=NT> idx_obs; 
  array[NT_mis] int<lower=1, upper=NT> idx_mis;

  // climate forcings
  vector[NT] G; // greenhouse gases
  vector[NT] S; // solar irradiance
  vector[NT] V; // volcanic activity

  vector[NT_prj] G_prj; // projected greenhouse gases (for future projections)
  vector[NT_prj] S_prj; // projected solar irradiance (for future projections)
  real vol_coef; // coefficient to scale the simulated volcanic activity for future projections (to match historical scale)

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
  vector[NT_obs] y_ins_fitted;  // fitted temperature values for the instrumental period (for posterior predictive checks)
  real sigma_y_ins;             // effective noise for the fitted values during the instrumental period (combining process and proxy noise)
  vector[NT_prj] V_prj = rep_vector(0, NT_prj); // future volcanic forcing (for projections)
  vector[NT_prj] y_prj;      // future temperature projections (for projections)

  // y mean (conditional on realized y[t-1]) 
  mu_y[1] = mu_forcing[1];
  for (t in 2:NT)
    mu_y[t] = phi_T * y[t-1] + mu_forcing[t];

  // y fitted values for the observed period (for PPCs)
  sigma_y_ins = alpha1^2/sigma_z^2;
  sigma_y_ins = sigma_y_ins + (1/sigma_y^2);
  sigma_y_ins = inv_sqrt(sigma_y_ins);

  for (t in 1:NT_obs){
    real mu_help;
    real mu_y_help;
    if (t == 1){
      mu_y_help = phi_T * y[NT_mis] + mu_forcing[idx_obs[t]];
      mu_help = (mu_y_help/sigma_y^2) + (alpha1 * (z[idx_obs[t]] - phi_R*z[NT_mis])/sigma_z^2); // CHECK THIS LINE
    }
    else {
      mu_y_help = phi_T * y_ins_fitted[t-1] + mu_forcing[idx_obs[t]];
      mu_help = (mu_y_help/sigma_y^2) + (alpha1 * (z[idx_obs[t]] - phi_R*z[idx_obs[t-1]])/sigma_z^2);
    }
    mu_help = mu_help * sigma_y_ins^2;
    y_ins_fitted[t] = normal_rng(mu_help, sigma_y_ins);
  }

  // v_simulate (simulate volcanism for y projections), model global volcanism (AOD) as a marked gamma process with a decay function.
  if (NT_prj > 0){
    vector[NT_prj] spikes_future = rep_vector(0, NT_prj);
    real baseV = 0.0; // baseline volcanic activity (can be adjusted)
    int current_spike_idx;
    int spike_0;
    int first_spike;
    int keep_going = 1;
    vector[8] Vdecay = [
      1.0,
      0.63871871,
      0.23100898,
      0.0786814,
      0.02663635,
      0.00928752,
      0.0031552,
      0.00209627
    ]';

    // NOTE: Per STAN documentation, gamma_rng returns a gamma random variate for the given shape 
    // and inverse scale (i.e., rate) parameters.

    spike_0 = -7;

    first_spike = spike_0 + to_int(gamma_rng(2.16, 0.2));
    while (first_spike < 1){
      first_spike = spike_0 + to_int(gamma_rng(2.16, 0.2));
    }

    current_spike_idx = first_spike;

    if (current_spike_idx <= NT_prj) {
      spikes_future[current_spike_idx] = gamma_rng(1.08, 19.5); // Amplitude of the spike
    }

    while (keep_going){
      int next_spike_time;
      next_spike_time = current_spike_idx + to_int(gamma_rng(2.16, 0.2)); // CHECK RATE VERSUS SCALE PARAMETERIZATION!!!
      if (next_spike_time < NT_prj){
        current_spike_idx = next_spike_time;
        spikes_future[current_spike_idx] = gamma_rng(1.08, 19.5); // CHECK RATE VERSUS SCALE PARAMETERIZATION!!!
      }
      else {
        keep_going = 0;
      }
    }

    V_prj = V_prj + rep_vector(baseV, NT_prj); // add baseline volcanic activity to all future years
    for (t in 1:NT_prj) {
      for (i in 0:7) {
        if (t - i >= 1) {
          V_prj[t] += spikes_future[t - i] * Vdecay[i + 1];
        }
      }
    }

    // normalize the simulated volcanic forcing to match the scale of historical volcanic forcing
    for (t in 1:NT_prj) {
      V_prj[t] = -abs(vol_coef) * (1 - exp(-V_prj[t]));
      }
  }

  // future temperature projections (for projections)
  if (NT_prj > 0){
    vector[NT_prj] mu_forcing_prj;
    mu_forcing_prj = beta0 + betaG * G_prj + betaS * S_prj + betaV * V_prj;

    // project future temperatures using the AR(1) structure
    y_prj[1] = normal_rng(phi_T * y[NT] + mu_forcing_prj[1], sigma_y);
    for (t in 2:NT_prj) {
      y_prj[t] = normal_rng(phi_T * y_prj[t-1] + mu_forcing_prj[t], sigma_y);
    }
  }

}
data {
  int                    n_priors;        // n priors
  matrix[n_priors,2]     priors;          // rows: priors; cols: p1, p2
  int                    n_years;         // n years. 24. SAME live, dead.
  int                    n_obs;          // number of spawner count obs
  
  int                    n_locations;    // number of spawning locations
  int<lower=1, upper=n_locations> location[n_obs]; //index for location
  int<lower=0, upper=1> is_valid[n_years, n_locations]; //1 if location is valid in that year 0 otherwise
  
  vector[n_locations] prior_weights;// dirichlet prior for relative abundance of fish in each location
  
  int<lower=1, upper=n_years>     year[n_obs];
  int                    day[n_obs];
  int                    live_counts[n_obs];
}
parameters {
  //run proportions for the different locations
  matrix[n_years, n_locations - 1] run_prop_offset;  // relative logits for locations 2..K
  vector[n_locations - 1] run_prop_mean;
  real<lower=0> run_prop_sigma;
  
  // group distributions (years)
  vector[n_locations] timing_mu_loc;         // average per location
  real<lower=0> timing_sigma;            // variation across years
  
  vector[n_locations]       spread_arrive_mu;     // stdev entry
  real<lower=0> spread_arrive_sigma;     
  
  vector[n_locations]       spread_death_mu;      // stdev death
  real<lower=0> spread_death_sigma;     
  
  // simple
  vector[n_locations]     residence_raw;  // lag entry to exit
  real<lower=0>         phi_live;   // dispersion parameter on live fish sampling
  
  // instance in group, multi-level
  vector [n_years] log_run;
  matrix[n_years, n_locations] timing_z;
  matrix[n_years, n_locations] spread_arrive_z;
  matrix[n_years, n_locations] spread_death_z;
}

transformed parameters{
  simplex[n_locations] run_prop[n_years];
  
  for (y in 1:n_years) {
    vector[n_locations] raw;
    raw[1] = is_valid[y, 1] == 1 ? 0 : -1e6;
    
    for (l in 2:n_locations) {
      if (is_valid[y, l] == 1) {
        raw[l] = run_prop_offset[y, l - 1];
      } else {
        raw[l] = -1e6;
      }
    }
    
    run_prop[y] = softmax(raw);
  }
  
  vector[n_locations] log_run_prop_masked[n_years];
  
  for (y in 1:n_years) {
    log_run_prop_masked[y] = rep_vector(log(1e-6), n_locations);
    for (l in 1:n_locations)
    if (is_valid[y, l] == 1)
    log_run_prop_masked[y, l] = log1p(run_prop[y, l] * 1e4) - log(1e4);
  }
  
  //non-centered priors
  matrix[n_years, n_locations] timing;
  matrix[n_years, n_locations] spread_arrive;
  matrix[n_years, n_locations] spread_death;
  
  
  for (y in 1:n_years){
    for (l in 1:n_locations){
      timing[y, l] = timing_mu_loc[l] + timing_sigma * timing_z[y, l];
      spread_arrive[y, l] = exp(spread_arrive_mu[l] + spread_arrive_sigma * spread_arrive_z[y,l]) + 1;
      spread_death[y, l] = exp(spread_death_mu[l] + spread_death_sigma * spread_death_z[y,l]) + 1;
    }
  }
  
  vector[n_locations] residence;
  residence = exp(residence_raw) + 1;
  
  //process model
  // live
  vector[n_obs] mu_live; // predictions from model
  for(i in 1:n_obs){
    real mean_arrival = timing[year[i], location[i]];
    real entered = normal_cdf(day[i], mean_arrival, spread_arrive[year[i], location[i]]);  
    real exited  = normal_cdf(day[i], mean_arrival + residence[location[i]], spread_death[year[i], location[i]]);
    mu_live[i] = exp(
      log_run[year[i]] +
      log_run_prop_masked[year[i], location[i]] +
      log1p(entered * (1 - exited) * 1e4) - log(1e4)
      );
  }
}

model {
  // local variables
  // group prior PDDs
  
  timing_mu_loc ~ normal(priors[3,1], priors[3,2]); 
  timing_sigma ~ exponential(1);                
  to_vector(timing_z) ~ normal(0, 1);             
  
  spread_arrive_mu  ~ normal( priors[6,1], priors[6,2]);  // 
  spread_arrive_sigma ~ exponential(priors[9,2]);
  to_vector(spread_arrive_z) ~ normal(0, 1);
  
  spread_death_mu   ~ normal( priors[6,1], priors[6,2]);  // 
  spread_death_sigma ~ exponential(priors[9,2]);
  to_vector(spread_death_z) ~ normal(0, 1);
  
  // simple PDD
  log_run ~ normal(priors[1,1], priors[1,2]);  // mean  run log
  residence_raw  ~ normal( priors[7,1], priors[7,2]);  // same all years
  log(phi_live)   ~ normal( priors[8,1], priors[8,2]);
  
  for (l in 1:(n_locations - 1))
  run_prop_mean[l] ~ normal(log(prior_weights[l + 1]) - log(prior_weights[1]), 0.5);
  
  run_prop_sigma ~ exponential(2);
  
  for (y in 1:n_years)
  for (l in 1:(n_locations - 1))
  run_prop_offset[y, l] ~ normal(run_prop_mean[l], run_prop_sigma);
  
  //likelihood
  target += neg_binomial_2_lpmf(live_counts | mu_live, phi_live);
  //target +=  normal_lpdf(log_dead_pred | dead_obs_day[,2], samp_err_dead);
}

generated quantities{
  vector[n_obs] log_lik;
  for(i in 1:n_obs)
  log_lik[i] = neg_binomial_2_lpmf(live_counts[i] | mu_live[i], phi_live);
}

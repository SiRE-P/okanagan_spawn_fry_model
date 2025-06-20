functions{
  vector merge_missing(array[] int miss_indexes, vector x_obs, vector x_miss) {
    int N = dims(x_obs)[1];
    int N_miss = dims(x_miss)[1];
    vector[N] merged;
    merged = x_obs;
    for(i in 1:N_miss)
    merged[ miss_indexes[i] ] = x_miss[i];
    return merged;
  }
}

data {
  int<lower=0> Y; //years
  int<lower=0> N; //number of observations
  int<lower=0> D; //number of days to be modelled
  int fry_obs[N]; //observed fry
  vector[N] soak_time; //soak time in minutes
  int year[N]; //year
  vector[N] hour; //hour set started - standardized 
  int day[N]; //observation day index
  vector[Y] new_moon_date; //date of the new moon, standardized
  vector[Y] ATU; //ATUs on day 100
  vector[N] dusk; //dusk on day of observation
  
  int<lower=1> inc_days;
  array[Y, inc_days]  real<lower=0> incubation_flow;
  real scour_threshold_prior_mu;
  real<lower=0> scour_threshold_prior_sigma;
  
  int freshet_days;
  array[Y, freshet_days]  real<lower=0> freshet_flow;
  real freshet_threshold_prior_mu;
  real<lower=0> freshet_threshold_prior_sigma;
  
  vector[Y] spawner_ln_est;
  vector[Y] spawner_ln_sd;
  
  vector[D] day_std; //sequence of all days to be predicted
  
  real sigma_sf_prior;
  real alpha_sf_prior;
  real beta_sf_prior;
  real alpha_sf_sigma_prior;
  real beta_sf_sigma_prior;
}

parameters {
  vector<lower=1>[Y] spawners;
  
  real beta0;
  real<lower=0> theta_sf;
  real alpha0; 
  real<lower=0> sigma_sf;
  vector[Y] total_fry_ln;
  
  real<lower=0> scour_threshold;  // estimated breakpoint
  real<lower=0> scour_transition_slope; 
  
  real<lower=0> freshet_threshold;  // estimated breakpoint
  real<lower=0> freshet_transition_slope; 
  
  real a_scour;
  real b_freshet;
  real sf_ATU;
  
  real<lower=0> soak_b;
  real b_moon;
  real b_ATU;
  real b_moonxATU;
  real b_dusk;
  
  real hour_peak_mu;
  real hour_sd_mu;
  real day_peak_mu;
  real day_sd_mu;
  
  real<lower=0> hour_peak_sigma;
  real<lower=0> hour_sd_sigma;
  real<lower=0> day_peak_sigma;
  real<lower=0> day_sd_sigma;
  
  vector[Y] hour_peak_z;
  vector[Y] hour_sd_z;
  vector[Y] day_peak_z;
  vector[Y] day_sd_z;
  
  real<lower=0> emerg_obs_error;
  
}
transformed parameters{
  
  //non-centered priors
  
  vector[Y] hour_sd;
  vector[Y] day_peak;
  vector[Y] day_sd;
  vector[Y] scour_count;
  vector[Y] fresh_count;
  
  hour_sd = exp(hour_sd_mu + hour_sd_sigma * hour_sd_z);
  day_peak = day_peak_mu + day_peak_sigma * day_peak_z + 
  b_moon * new_moon_date + 
  b_ATU * ATU +
  b_moonxATU * (new_moon_date .* ATU);
  day_sd = exp(day_sd_mu + day_sd_sigma * day_sd_z);
  
  
  
  //year specific alphas
  vector[Y] alpha_sf;
  vector[Y] beta_sf;
  for (y in 1:Y) {
        scour_count[y] = 0;
        fresh_count[y] = 0;
    for (d in 1:inc_days) {
      scour_count[y] += inv_logit(scour_transition_slope * (incubation_flow[y, d] - scour_threshold));    
    }
    for (d in 1:freshet_days) {
      fresh_count[y] += inv_logit(freshet_transition_slope * (freshet_flow[y, d] - freshet_threshold));    
    }
    
    alpha_sf[y] = exp(alpha0
    + a_scour * scour_count[y]/5 +
    + sf_ATU * ATU[y]);
    
    beta_sf[y] = exp(beta0 + b_freshet * fresh_count[y]/30);
  }
  
  vector[Y] peak_fry;
  for (y in 1:Y){
    peak_fry[y] = total_fry_ln[y] - log((sqrt(2 * pi()) * day_sd[y] * sqrt(2 * pi()) * hour_sd[y]));
  }
  
  matrix[Y, D] emerging_fry_ln;
  for (y in 1:Y){
    for (d in 1:D){
      emerging_fry_ln[y, d] = peak_fry[y] - ((day_std[d] - day_peak[y])^2) / (2 * day_sd[y]^2);
    }
  }
}

model {
  //spawner fry beverton-holt
  
  alpha0 ~ normal(alpha_sf_prior, alpha_sf_sigma_prior);
  beta0 ~ normal(beta_sf_prior, beta_sf_sigma_prior);
  theta_sf ~ normal(1, 0.1);
  
  scour_threshold ~ normal(scour_threshold_prior_mu, scour_threshold_prior_sigma);
  freshet_threshold ~ normal(freshet_threshold_prior_mu, freshet_threshold_prior_sigma);
  scour_transition_slope ~ lognormal(1.5, 0.6); // median ~4.5, 95% range ≈ [1.1, 18]
  freshet_transition_slope ~ lognormal(1.5, 0.6); // median ~4.5, 95% range ≈ [1.1, 18]
  a_scour ~ normal(0, 0.5);
  b_freshet ~ normal(0, 0.5);
  sf_ATU ~ normal(0, 0.5);
  
  spawners ~ lognormal(spawner_ln_est, spawner_ln_sd);
  
  sigma_sf ~ exponential(sigma_sf_prior);
  for (y in 1:Y){
    real fry_mu_ln = log((alpha_sf[y] * spawners[y])/(1 + (beta_sf[y] * spawners[y]/1e05)^theta_sf));
    total_fry_ln[y] ~ normal(fry_mu_ln, sigma_sf);
  }
  
  //fry observation model
  soak_b ~ exponential(1);
  b_moon ~ normal(0, 0.5);
  b_ATU ~ normal(0, 0.5);
  b_moonxATU ~ normal(0, 0.5);
  b_dusk ~ normal(0, 0.5);
  
  hour_peak_mu ~ normal(0, 0.5);
  hour_sd_mu ~ normal(-0.5, 0.5);
  hour_peak_sigma ~ exponential(10);
  hour_sd_sigma ~ exponential(10);
  hour_peak_z ~ normal(0, 1);
  hour_sd_z ~ normal(0, 1);
  
  day_peak_mu ~ normal(0, 0.5);
  day_sd_mu ~ normal(-0.5, 0.5);
  day_peak_sigma ~ exponential(1);
  day_sd_sigma ~ exponential(1);
  day_peak_z ~ normal(0, 1);
  day_sd_z ~ normal(0, 1);
  
  emerg_obs_error ~ exponential(1);
  
  for(i in 1:N){
    real hour_peak = hour_peak_mu + hour_peak_sigma * hour_peak_z[year[i]] +
    b_dusk * dusk[i];
    real obs_offset = -((hour[i] - hour_peak)^2) / (2*hour_sd[year[i]]^2) + soak_b * soak_time[i];
    real emerge_mu = exp(emerging_fry_ln[year[i], day[i]] + obs_offset);
    fry_obs[i] ~ neg_binomial_2(emerge_mu, 1/emerg_obs_error);
  }
}
generated quantities{
  vector[Y] log_lik;
  for (y in 1:Y) {
    real fry_mu_ln = log((alpha_sf[y] * spawners[y])/(1 + (beta_sf[y] * spawners[y]/1e5)^theta_sf));
    log_lik[y] = normal_lpdf(total_fry_ln[y] | fry_mu_ln, sigma_sf);
  }
}

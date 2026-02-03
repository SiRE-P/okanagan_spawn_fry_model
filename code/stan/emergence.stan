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
  array[N] int fry_obs; //observed fry
  
  int N_missing_volume;
  array[N_missing_volume] int volume_missidx; 
  vector[N] volume;
  vector[N] sample_flow;
  
  vector[N] soak_time; //soak time in minutes
  array[N] int year; //year
  vector[N] hour; //hour set started - standardized 
  array[N] int day; //observation day index
  vector[Y] new_moon_date; //date of the new moon, standardized
  vector[N] dusk; //dusk on day of observation
  array[Y] int exceeded_FWMT; //index for whether or not FWMT range was exceeded 1 for no, 2 for yes
  
  vector[Y] spawner_ln_est;
  vector[Y] spawner_ln_sd;
  
  int N_missing_ATU;
  array[N_missing_ATU] int ATU_missidx;
  vector[Y] ATU; //ATUs on day 100
  
  int N_missing_thermal;
  array[N_missing_thermal] int thermal_missidx;
  vector[Y] thermal_onset;
  vector[Y] thermal_dur_resid;
  
  vector[Y] malott_june_air;
  
  vector[D] day_std; //sequence of all days to be predicted
  
  int freshet_days;
  array[Y, freshet_days]  real<lower=0> freshet_flow;
  real freshet_transition_slope;
  real freshet_threshold_prior_mu;
  real<lower=0> freshet_threshold_prior_sigma;

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
  
  real a_FWMT;
  real sf_ATU;
  
  real<lower=0> freshet_threshold; 

  real a_effective_spawners;
  real b_therm_onset_raw;
  real b_therm_dur_raw;
  real<lower=0> thermal_barrier_transition_width_sd; 
  
  real a_onset;
  real b_air_onset;
  real<lower=0> onset_sigma;
  
  real b_freshet;
  
  real a_volume;
  real soak_b;
  real b_sample_flow;
  real<lower=0> sigma_volume;
  
  real b_volume;
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
  
  real<lower=1> emerg_phi;
  
  //missing data
  vector<lower=0>[N_missing_volume] volume_impute;
  vector[N_missing_ATU] ATU_impute;
  vector[N_missing_thermal] thermal_onset_impute;
  vector[N_missing_thermal] thermal_dur_resid_impute;
}
transformed parameters{
  //missing variables
  vector[N] volume_merge;
  volume_merge = merge_missing(volume_missidx, to_vector(volume), volume_impute);

  vector[Y] ATU_merge;
  ATU_merge = merge_missing(ATU_missidx, to_vector(ATU), ATU_impute);

  vector[Y] thermal_onset_merge;
  thermal_onset_merge = merge_missing(thermal_missidx, to_vector(thermal_onset), thermal_onset_impute);

  vector[Y] thermal_dur_resid_merge;
  thermal_dur_resid_merge = merge_missing(thermal_missidx, to_vector(thermal_dur_resid), thermal_dur_resid_impute);
  
  
  //non-centered priors
  vector[Y] hour_sd;
  vector[Y] day_peak;
  vector[Y] day_sd;
  vector[Y] fresh_days;
  vector[Y] fresh_days_scaled;
  vector[Y] effective_spawners;
  vector[Y] eff_spawner_prop;
  
  hour_sd = exp(hour_sd_mu + hour_sd_sigma * hour_sd_z);
  day_peak = day_peak_mu + day_peak_sigma * day_peak_z + 
  b_moon * new_moon_date + 
  b_ATU * ATU_merge +
  b_moonxATU * (new_moon_date .* ATU_merge);
  day_sd = exp(day_sd_mu + day_sd_sigma * day_sd_z);
  
  //effective spawner paramaters scaled to SD units of onset date
  real slope_scale = 5.9 / thermal_barrier_transition_width_sd; //5.9 is the logit scale distance between 5% and 95%
  real b_therm_onset = slope_scale * b_therm_onset_raw;
  real b_therm_dur   = slope_scale * b_therm_dur_raw;
  
  //year specific alphas
  vector[Y] alpha_sf;
  vector[Y] beta_sf;
  for (y in 1:Y) {
  fresh_days[y] = 0;
    for (d in 1:freshet_days) {
      fresh_days[y] += inv_logit(freshet_transition_slope * (freshet_flow[y, d] - freshet_threshold));
    }    
    
    eff_spawner_prop[y] = inv_logit(a_effective_spawners + b_therm_onset * thermal_onset_merge[y] + b_therm_dur * thermal_dur_resid_merge[y]);
    effective_spawners[y] = eff_spawner_prop[y] * spawners[y];
    
    alpha_sf[y] = exp(alpha0
    + a_FWMT * exceeded_FWMT[y]
    + sf_ATU * ATU[y]
    );
  }
    fresh_days_scaled = fresh_days/sd(fresh_days);
    beta_sf = exp(beta0 + b_freshet * fresh_days_scaled); // 57 = sd(rowSums(freshet_mat>28.5))

  
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
  ATU_impute ~ std_normal();
  
  //spawner fry beverton-holt
  
  alpha0 ~ normal(alpha_sf_prior, alpha_sf_sigma_prior);
  beta0 ~ normal(beta_sf_prior, beta_sf_sigma_prior);
  theta_sf ~ normal(1, 0.2);
  
  a_FWMT ~ normal(0, 0.5);
  sf_ATU ~ normal(0, 0.5);

  freshet_threshold ~ normal(freshet_threshold_prior_mu, freshet_threshold_prior_sigma);

  a_effective_spawners ~ normal(0, 0.25); 
  thermal_barrier_transition_width_sd ~ lognormal(log(2.5), 0.4);
  b_therm_dur_raw ~ normal(0, 0.5);
  b_therm_onset_raw ~ normal(0, 0.5);
  b_freshet ~ normal(0, 0.5);
  
  spawners ~ lognormal(spawner_ln_est, spawner_ln_sd);
  
  thermal_dur_resid_impute ~  std_normal();
  //thermal_onset_impute ~ std_normal();
  a_onset ~ normal(0, 0.1);
  b_air_onset ~ normal(-0.5, 0.5);
  onset_sigma ~ exponential(1);
  thermal_onset_merge ~ normal(a_onset + b_air_onset * malott_june_air, onset_sigma);
  
  
  a_volume ~ normal(0,2);
  soak_b ~ normal(1, 0.5);
  b_sample_flow ~ normal(1, 0.5);
  sigma_volume ~ exponential(10);
  volume_merge ~ lognormal(a_volume + soak_b * log(soak_time) + b_sample_flow * log(sample_flow), sigma_volume);

  
  sigma_sf ~ exponential(sigma_sf_prior);
  for (y in 1:Y){
    real fry_mu_ln = log((alpha_sf[y] * effective_spawners[y])/(1 + (beta_sf[y] * spawners[y]/1e5)^theta_sf));
    total_fry_ln[y] ~ normal(fry_mu_ln, sigma_sf);
  }
  
  //fry observation model
  b_volume ~ normal(1, 0.5);
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
  
  emerg_phi ~ lognormal(log(20),0.5);
  
  for(i in 1:N){
    real hour_peak = hour_peak_mu + hour_peak_sigma * hour_peak_z[year[i]] +
    b_dusk * dusk[i];
    real obs_offset = -((hour[i] - hour_peak)^2) / (2*hour_sd[year[i]]^2) + b_volume * (log(volume_merge[i]) - log(100));
    real emerge_mu = exp(emerging_fry_ln[year[i], day[i]] + obs_offset);
    fry_obs[i] ~ neg_binomial_2(emerge_mu, emerg_phi);
  }
}

generated quantities{
  vector[N] log_lik;
  for(i in 1:N){
    real hour_peak = hour_peak_mu + hour_peak_sigma * hour_peak_z[year[i]] + b_dusk * dusk[i];
    real obs_offset = -((hour[i] - hour_peak)^2) / (2*hour_sd[year[i]]^2) + b_volume * (log(volume_merge[i]) - log(100));
    real emerge_mu = exp(emerging_fry_ln[year[i], day[i]] + obs_offset);
    log_lik[i] = neg_binomial_2_lpmf(fry_obs[i] | emerge_mu, emerg_phi);
  }
}

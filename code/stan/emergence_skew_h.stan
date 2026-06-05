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
  //vector[Y] new_moon_date; //date of the new moon, standardized
  //vector[N] dusk; //dusk on day of observation
  array[Y] int exceeded_FWMT; //index for whether or not FWMT range was exceeded 0 for no, 1 for yes
  
  vector[Y] spawner_ln_est;
  vector[Y] spawner_ln_sd;
  
  int N_missing_ATU;
  array[N_missing_ATU] int ATU_missidx;
  vector[Y] ATU; //ATUs on day 100
  
  int N_missing_thermal;
  array[N_missing_thermal] int thermal_missidx;
  vector[Y] thermal_onset;
  real onset_mean;
  real onset_sd;
  real onset_low_bound;
  vector[Y] thermal_dur_resid;
  
  vector[Y] malott_june_air;
  
  vector[D] day_std; //sequence of all days to be predicted
  
  int freshet_days;
  array[Y, freshet_days]  real<lower=0> freshet_flow;
  real freshet_transition_slope;
  real freshet_threshold;
  
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
  
  real hour_peak;
  real<lower=0> hour_sd;
  
  vector[3] emerge_phenology_mu;
  vector<lower=0>[3] emerge_phenology_sigma;
  cholesky_factor_corr[3] emerge_phenology_L_corr;
  matrix[3, Y] emerge_phenology_z;
  
  real b_ATU_peak; //ATU effect on emergence peak
  
  real<lower=1> emerg_phi;
  
  //missing data
  vector<lower=0>[N_missing_volume] volume_impute;
  vector[N_missing_ATU] ATU_impute;
  vector<lower=onset_low_bound>[N_missing_thermal] thermal_onset_impute;
  vector[N_missing_thermal] thermal_dur_resid_impute;
  
  //hyperparamters for imputation of duration  
  real xi_dur;
  real<lower=0> omega_dur;
  real alpha_dur;
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
  
  //standardized thermal onset
  vector[Y] thermal_onset_scaled;
  thermal_onset_scaled = (thermal_onset_merge - onset_mean) / onset_sd;
  
  //non-centered priors
  matrix[Y, 3] emerge_phenology;
  vector[Y] day_peak;
  vector[Y] day_sd;
  vector[Y] day_skew;
  
  for (y in 1:Y) {
    vector[3] mu_y;
    mu_y[1] = emerge_phenology_mu[1]+ b_ATU_peak * ATU_merge[y];
    mu_y[2] = emerge_phenology_mu[2];
    mu_y[3] = emerge_phenology_mu[3];
    
    emerge_phenology[y] = to_row_vector(mu_y + 
    diag_pre_multiply(emerge_phenology_sigma, emerge_phenology_L_corr) * emerge_phenology_z[, y]);
    
    day_peak[y] = emerge_phenology[y, 1];
    day_sd[y]   = exp(emerge_phenology[y, 2]);
    day_skew[y] = emerge_phenology[y, 3];
  }
  
  vector[Y] fresh_days;
  vector[Y] fresh_days_scaled;
  vector[Y] effective_spawners;
  vector[Y] eff_spawner_prop;
  
  
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
    
    eff_spawner_prop[y] = inv_logit(a_effective_spawners + b_therm_onset * thermal_onset_scaled[y] + b_therm_dur * thermal_dur_resid_merge[y]);
    effective_spawners[y] = eff_spawner_prop[y] * spawners[y];
    
    alpha_sf[y] = exp(alpha0
    + a_FWMT * exceeded_FWMT[y]
    + sf_ATU * ATU_merge[y]
    );
  }
  fresh_days_scaled = fresh_days/sd(fresh_days);
  beta_sf = exp(beta0 + b_freshet * fresh_days_scaled); // 57 = sd(rowSums(freshet_mat>28.5))
  
  
  matrix[Y, D] day_kernel_ln;
  vector[Y] day_kernel_log_norm;
  
  for (y in 1:Y) {
    for (d in 1:D) {
      real z = (day_std[d] - day_peak[y]) / day_sd[y];
      day_kernel_ln[y, d] = -0.5 * square(z) + normal_lcdf(day_skew[y] * z | 0, 1);
    }
    day_kernel_log_norm[y] = log_sum_exp(day_kernel_ln[y]);
  }
  
  matrix[Y, D] emerging_fry_ln;
  for (y in 1:Y){
    for (d in 1:D){
      emerging_fry_ln[y, d] = total_fry_ln[y] + day_kernel_ln[y, d] - day_kernel_log_norm[y];
    }
  }
}

model {
  ATU_impute ~ std_normal();
  
  //spawner fry beverton-holt
  
  alpha0 ~ normal(alpha_sf_prior, alpha_sf_sigma_prior);
  beta0 ~ normal(beta_sf_prior, beta_sf_sigma_prior);
  theta_sf ~ normal(1, 0.1);
  
  a_FWMT ~ normal(0, 0.5);
  sf_ATU ~ normal(0, 0.5);
  
  a_effective_spawners ~ normal(0, 0.25); 
  thermal_barrier_transition_width_sd ~ lognormal(log(2.5), 0.5);
  b_therm_dur_raw ~ normal(0, 0.5);
  b_therm_onset_raw ~ normal(0, 0.5);
  b_freshet ~ normal(0, 0.5);
  
  spawners ~ lognormal(spawner_ln_est, spawner_ln_sd);
  
  // priors informed by observed standardized data
  xi_dur    ~ normal(0.5, 0.4);    
  omega_dur ~ normal(1.30, 0.2); 
  alpha_dur ~ normal(-2, 1);
  
  thermal_dur_resid_merge ~ skew_normal(xi_dur, omega_dur, alpha_dur);
  
  a_onset ~ normal(0, 0.1);
  b_air_onset ~ normal(-0.5, 0.5);
  onset_sigma ~ exponential(1);
  thermal_onset_scaled ~ normal(a_onset + b_air_onset * malott_june_air, onset_sigma);
  
  
  a_volume ~ normal(0,1);
  soak_b ~ normal(1, 0.5);
  b_sample_flow ~ normal(1, 0.5);
  sigma_volume ~ lognormal(log(0.3), 0.4);
  volume_merge ~ lognormal(a_volume + soak_b * log(soak_time) + b_sample_flow * log(sample_flow), sigma_volume);
  
  
  sigma_sf ~ exponential(sigma_sf_prior);
  for (y in 1:Y){
    real fry_mu_ln = log((alpha_sf[y] * effective_spawners[y])/(1 + (beta_sf[y] * spawners[y]/1e5)^theta_sf));
    total_fry_ln[y] ~ normal(fry_mu_ln, sigma_sf);
  }
  
  //fry observation model
  b_volume ~ normal(1, 0.5);
  //b_dusk ~ normal(0, 0.5);
  
  hour_peak ~ normal(-0.5, 0.5);
  hour_sd ~ lognormal(log(1), 0.4);
  
  
  // Emergence phenology
  emerge_phenology_mu[1] ~ normal(0, 0.5);      
  emerge_phenology_mu[2] ~ normal(-0.5, 0.5);
  emerge_phenology_mu[3] ~ normal(-2, 1);
  
  emerge_phenology_sigma[1] ~ lognormal(log(0.3), 0.4);
  emerge_phenology_sigma[2] ~ lognormal(log(0.3), 0.4);
  emerge_phenology_sigma[3] ~ lognormal(log(2), 0.4);
  
  emerge_phenology_L_corr ~ lkj_corr_cholesky(2);
  
  to_vector(emerge_phenology_z) ~ normal(0, 1);
  
  b_ATU_peak   ~ normal(-0.2, 0.5);
  
  emerg_phi ~ lognormal(log(5),0.5);
  
  for(i in 1:N){
    real obs_offset = -((hour[i] - hour_peak)^2) / (2*hour_sd^2) + b_volume * (log(volume_merge[i]) - log(100));
    real emerge_mu = exp(emerging_fry_ln[year[i], day[i]] + obs_offset);
    fry_obs[i] ~ neg_binomial_2(emerge_mu, emerg_phi);
  }
}

generated quantities{
  vector[N] log_lik;
  for(i in 1:N){
    real obs_offset = -((hour[i] - hour_peak)^2) / (2*hour_sd^2) + b_volume * (log(volume_merge[i]) - log(100));
    real emerge_mu = exp(emerging_fry_ln[year[i], day[i]] + obs_offset);
    log_lik[i] = neg_binomial_2_lpmf(fry_obs[i] | emerge_mu, emerg_phi);
  }
}

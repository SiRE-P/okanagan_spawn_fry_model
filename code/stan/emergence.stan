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
  int exceeded_FWMT[Y]; //index for whether or not FWMT range was exceeded 1 for no, 2 for yes

  int N_missing_spawners;
  array[N_missing_spawners]  int spawners_missidx;
  vector[Y] spawners;
  real log_spawners_mean;
  real log_spawners_sd;
  
  vector[D] day_std; //sequence of all days to be predicted
  
  real sigma_sf_prior;
  real beta_sf_prior;
  real eq_ratio_prior;
  real beta_sf_sigma_prior;
  real eq_ratio_sigma_prior;
}

parameters {
  real<lower=0> beta_sf;
  real<lower=0> eq_ratio; 
  real<lower=0> sigma_sf;
  vector[Y] total_fry_ln;
  
  real a_FWMT;
  real sf_ATU;
  
  real<lower=0> soak_b;
  real b_moon;
  real b_ATU;
  real b_moonxATU;
  
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
  
  //missing data
  vector<lower=0>[N_missing_spawners] spawners_impute;
}
transformed parameters{
  //missing variables
  vector[Y] spawners_merge;
  spawners_merge = merge_missing(spawners_missidx, to_vector(spawners), spawners_impute);
  
  //non-centered priors
  real alpha0 = log(eq_ratio * beta_sf);

  vector[Y] hour_peak;
  vector[Y] hour_sd;
  vector[Y] day_peak;
  vector[Y] day_sd;
  
  hour_peak = hour_peak_mu + hour_peak_sigma * hour_peak_z;
  hour_sd = exp(hour_sd_mu + hour_sd_sigma * hour_sd_z);
  day_peak = day_peak_mu + day_peak_sigma * day_peak_z + 
        b_moon * new_moon_date + 
        b_ATU * ATU +
        b_moonxATU * (new_moon_date .* ATU);
  day_sd = exp(day_sd_mu + day_sd_sigma * day_sd_z);
  
  //year specific alphas
  vector[Y] alpha_sf;
  for (y in 1:Y) {
    alpha_sf[y] = exp(alpha0
                      + a_FWMT * exceeded_FWMT[y]
                      + sf_ATU * ATU[y]);
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
  beta_sf ~ lognormal(beta_sf_prior, beta_sf_sigma_prior);
  eq_ratio ~ lognormal(eq_ratio_prior, eq_ratio_sigma_prior);
  
  a_FWMT ~ normal(0, 0.5);
  sf_ATU ~ normal(0, 0.5);

  spawners_merge ~ lognormal(log_spawners_mean, log_spawners_sd);

  sigma_sf ~ exponential(sigma_sf_prior);
  for (y in 1:Y){
    real fry_mu_ln = log((alpha_sf[y] * spawners_merge[y])/(1 + beta_sf * spawners_merge[y]/1e05));
    total_fry_ln[y] ~ normal(fry_mu_ln, sigma_sf);
  }
  
  //fry observation model
  soak_b ~ exponential(1);
  b_moon ~ normal(0, 0.5);
  b_ATU ~ normal(0, 0.5);
  b_moonxATU ~ normal(0, 0.5);
  
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
    real obs_offset = -((hour[i] - hour_peak[year[i]])^2) / (2*hour_sd[year[i]]^2) + soak_b * soak_time[i];
    real emerge_mu = exp(emerging_fry_ln[year[i], day[i]] + obs_offset);
    fry_obs[i] ~ neg_binomial_2(emerge_mu, 1/emerg_obs_error);
  }
}


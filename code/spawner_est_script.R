library(tidyverse)
library(tidybayes)
library(rstan)
library(loo)
library(patchwork)
library(ggrepel)
library(janitor)
library(rvest)

spawners <- read.csv("../../okanagan_data/2025-04-07 Draft/Spawn Timing Data/spawn_timing.csv") %>% 
  mutate(date = as.Date(date), year = year(date), yday = yday(date))

oso_sk <- spawners %>% 
  mutate(location = ifelse(location == "Channel", "VDS", location)) %>% 
  group_by(location, year, yday, date) %>% 
  filter(location != "Whole System") %>% 
  filter(species == "sk") %>% 
  summarise(
    live = if (all(is.na(live))) NA_real_ else sum(live, na.rm = TRUE),
    dead = if (all(is.na(dead))) NA_real_ else sum(dead, na.rm = TRUE)) %>% 
  arrange(date) %>% 
  ungroup() %>% 
  mutate(day_ind = as.numeric(factor(yday)), 
         year_ind = as.numeric(factor(year))) %>% 
  filter(!(is.na(live) & is.na(dead))) %>% 
  mutate(location = factor(location, levels = c("Index", "VDS", "Above McIntyre Dam"), ordered = TRUE)) %>% 
  #filter(!(location == "Above McIntyre Dam" & year == 2010)) %>% 
  filter(live>0)

oso_sk %>% 
  ggplot(aes(x = yday, y = live, color = location))+
  geom_point()+
  facet_wrap(~year)+
  scale_y_log10()

oso_sk %>% 
  ggplot(aes(x = yday, y = live, color = location))+
  geom_point()+
  facet_wrap(~year)

valid_id <- oso_sk %>% 
  group_by(year, location) %>% 
  summarise(valid = sum(live)>0) %>% 
  spread(key = location, value = valid, fill = 0) %>% 
  mutate(VDS = 1) %>% 
  ungroup() %>% 
  select(-year) %>% 
  data.matrix()

mean_prop_abund <- oso_sk %>% 
  group_by(year, location) %>% 
  summarise(live = sum(live)) %>% 
  group_by(year) %>% 
  mutate(live_prop = live/sum(live)) %>% 
  ungroup() %>% 
  group_by(location) %>% 
  summarise(live_prop = mean(live_prop))

#Wells counts####
Wells_data <- read_html("https://www.cbr.washington.edu/dart/cs/php/rpt/adult_annual.php?sc=1&outputFormat=html&proj=WEL&startdate=1%2F1&enddate=12%2F31&run=") %>% 
  html_table(fill = TRUE) %>% 
  .[[1]] %>% 
  filter(Project!="Project") %>% 
  clean_names() %>% 
  type.convert(as.is = TRUE) %>% 
  mutate(count_hrs = ifelse(year<1998, 16, 24), guess_16hrs = FALSE)

Wells16_24 <- readxl::read_excel("../../Ok_sock_lifecycle_model/data/24_hr_dam_counts/Difference between 16- & 24-hour Wells Sockeye counts (1998-2022).xlsx", skip = 2) %>% 
  clean_names()

Wells_data <- bind_rows(Wells_data, Wells16_24 %>% select(year, sockeye = x16_hour) %>% mutate(count_hrs = 16, project = "Wells", guess_16hrs = FALSE))


ggplot(Wells_data, aes(x = year, y = sockeye, color = factor(count_hrs)))+
  geom_point()

#sigmas on timing exp(20) - expect ~ 10 days variation 
#and on spread exp(40) - expect ~ 5 days variation

location_levels <- oso_sk %>% 
  select(location) %>% 
  unique() %>% 
  mutate(loc_ind = 1:n())

priors <- data.frame(prior = c("log_runs_mu", "log_runs_sigma", "timings_mu", "timings_sigma", "dead_frac", "spread", "residence_mu", "residuals", "residence_sigma", "sigma_r", "rho", "sigma_year"),
                     v1 = c(10,0,280-240,0,0,log(8-1),log(11-1), 1, 0, 0, 5, 0), 
                     v2 = c(1,10,5,1,2,0.2,0.3, 0.2, 40, 0.1, 5, 0.1))


sp_dat = list(n_priors = nrow(priors), 
              priors = data.matrix(priors[,-1]),
              n_years = max(oso_sk$year_ind),
              year = as.numeric(factor(oso_sk$year)),
              year_start = which(!duplicated(oso_sk$year)), 
              year_end = which(!duplicated(oso_sk$year, fromLast = TRUE)),
              n_locations = length(unique(oso_sk$location)),
              location = as.numeric(oso_sk$location),
              is_valid = valid_id,
              prior_weights = mean_prop_abund$live_prop,
              day = oso_sk$yday-240,
              n_obs = nrow(oso_sk),
              live_counts = oso_sk$live,
              dead_counts = oso_sk$dead)

m <- stan_model(file = "./code/stan/spawners.stan")
fit <- sampling(m, data = sp_dat, chains = 4, cores = 4)

worst_Rhat <- summary(fit)$summary %>% 
  as.data.frame() %>% 
  mutate(Rhat = round(Rhat, 3)) %>% 
  arrange(desc(Rhat))

worst_Rhat %>% 
  filter(n_eff>3) %>% 
  ggplot(aes(x = n_eff, y = Rhat))+
  geom_point()+
  geom_hline(yintercept = 1.01, lty = 2)+
  geom_vline(xintercept = 400, lty = 2)

traceplot(fit, pars = rownames(worst_Rhat)[1:20])
traceplot(fit, pars = c("timing_sigma", "timing_mu_loc", "spread_arrive_mu", "spread_arrive_sigma", "spread_death_mu", "spread_death_sigma"))

post <- extract(fit)

spread_draws(fit, residence[location]) %>% 
  ggplot(aes(x = location, y = residence))+
  geom_hline(yintercept = 11, lty = 2)+
  stat_pointinterval()

spread_draws(fit, run_prop[year, loc_ind]) %>% 
  left_join(location_levels) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  ggplot(aes(x = year, y = run_prop, color = factor(location), group = location))+
  stat_pointinterval()

spread_draws(fit, run_prop[year, loc_ind]) %>% 
  left_join(location_levels) %>% 
  left_join(spread_draws(fit, log_run[year])) %>% 
  mutate(run = exp(log_run) * run_prop) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  ggplot(aes(x = year, y = run))+
  stat_pointinterval()+
  facet_wrap(~location, scales = "free_y", ncol = 1)


#get summary stats for index sites for lining up with emergence data
spawner_est <- spread_draws(fit, run_prop[year, loc_ind]) %>% 
  left_join(location_levels) %>% 
  left_join(spread_draws(fit, log_run[year])) %>% 
  mutate(run = exp(log_run) * run_prop) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  filter(location == "Index") %>% 
  group_by(year) %>% 
  summarise(spawners = median(run), spawner_sd = sd(run))

write.csv(spawner_est, file = "./outputs/spawners_est.csv")

spread_draws(fit, log_run[year]) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  left_join(Wells_data %>% filter(count_hrs == 24) %>% select(year, Wells_count = sockeye)) %>% 
  ggplot(aes(x = Wells_count, y = exp(log_run), color = year>2010))+
  stat_pointinterval()+
  scale_color_discrete(name = "Skaha\npassage\nopen") +
  geom_text_repel(data = spread_draws(fit, log_run[year]) %>% 
                    mutate(year = year + min(spawners$year)-1) %>% 
                    left_join(Wells_data %>% filter(count_hrs == 24) %>% select(year, Wells_count = sockeye)) %>% 
                    group_by(year) %>% 
                    summarise(log_run = mean(log_run), Wells_count = mean(Wells_count)),
                  aes(label = year), color = 1)+
  geom_abline(slope = 1, intercept = 0, lty = 2)

spread_draws(fit, log_run[year]) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  ggplot(aes(x = year, y = exp(log_run)))+
  stat_pointinterval()+
  geom_point(data = spawners_old %>% filter(brood_year > 1999), aes(x = brood_year, y = spawners), color = "red")+
  
  spread_draws(fit, timing[year, loc_ind]) %>% 
  left_join(location_levels) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  filter(!(loc_ind == 3 & year < 2010)) %>% 
  ggplot(aes(x = year, y = timing, color = location, group = location))+
  stat_pointinterval()+
  
  spread_draws(fit, spread_arrive[year, loc_ind]) %>% 
  left_join(location_levels) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  filter(!(loc_ind == 3 & year < 2010)) %>% 
  ggplot(aes(x = year, y = spread_arrive, color = location, group = location))+
  stat_pointinterval()+
  
  spread_draws(fit, spread_death[year, loc_ind]) %>% 
  left_join(location_levels) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  filter(!(loc_ind == 3 & year < 2010)) %>% 
  ggplot(aes(x = year, y = spread_death, color = location, group = location))+
  stat_pointinterval()+
  plot_layout(ncol = 1)

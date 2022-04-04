## building a BYM route-level trend model for the BBS
setwd("C:/GitHub/BBS_iCAR_route_trends")
library(bbsBayes)
library(tidyverse)
library(cmdstanr)
# library(rstan)
# rstan_options(auto_write = TRUE, javascript = FALSE)
# library(shinystan)
library(sf)
library(spdep)
# library(doParallel)
# library(foreach)
 library(ggforce)
#library(tidybayes)
#source("functions/mungeCARdata4stan.R")
source("functions/neighbours_define.R") ## function to define neighbourhood relationships
source("functions/prepare-jags-data-alt.R") ## small alteration of the bbsBayes function
source("functions/get_basemap_function.R") ## loads one of the bbsBayes strata maps
source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output
## changes captured in a commit on Nov 20, 2020


# load and stratify CASW data ---------------------------------------------
#species = "Pacific Wren"
#species = "Barred Owl"
strat = "bbs_usgs"
model = "slope"


firstYear = 2004
lastYear = 2019

scope = "RangeWide"


species = "Blue-headed Vireo"
#species = "Dickcissel"

species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)


# SPECIES LOOP ------------------------------------------------------------

output_dir <- "output"


sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")

    load(sp_file) 
output_dir <- "output"
spp = "Stan_slope"

out_base <- paste0(species_f,spp,firstYear,"_",lastYear)

stan_data[["N_edges"]] = NULL
stan_data[["node1"]] = NULL
stan_data[["node2"]] = NULL

stan_data[["nstrata"]] = jags_data$nstrata
stan_data[["strat"]] = jags_data$strat



mod.file = "models/slope_NONiCAR_strata.stan"

## compile model
slope_model <- cmdstan_model(mod.file)

init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts,0,0.1),
                             alpha_raw = rnorm(stan_data$nroutes,0,0.1),
                             ALPHA = 0,
                             BETA = 0,
                             eta = 0,
                             obs_raw = rnorm(stan_data$nobservers,0,0.1),
                             sdnoise = 0.2,
                             sdobs = 0.1,
                             sdbeta_rand = runif(1,0.01,0.1),
                             beta_raw_rand = rnorm(stan_data$nstrata,0,0.01),
                             sdstrata = runif(1,0.01,0.1),
                             strata_raw = rnorm(stan_data$nstrata,0,0.01))}


slope_model <- cmdstan_model(mod.file)

slope_stanfit <- slope_model$sample(
  data=stan_data,
  refresh=200,
  chains=3, iter_sampling=1000,
  iter_warmup=1000,
  parallel_chains = 3,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 14,
  seed = 123,
  init = init_def,
  output_dir = output_dir,
  output_basename = out_base)

source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output

library(posterior)
stanf_df_stratum <- slope_stanfit$draws(format = "df")


conv_summ <- summarise_draws(stanf_df_stratum) %>% 
  mutate(species = species,
         model = out_base)


# extract trends and abundances -------------------------------------------

strat_df <- data.frame(route = jags_data$route,
                        strat = jags_data$strat,
                        stratum = jags_data$strat_name,
                        year = jags_data$r_year,
                        obs = jags_data$obser) %>% 
  group_by(strat,stratum,obs,route) %>%
  summarise(n_yr_obs_r = n(),
            .groups = "drop") %>% 
  group_by(strat,stratum,route) %>% 
  summarise(n_rt_obs = n(),
            mean_y_rt_obs = mean(n_yr_obs_r),
            max_y_rt_obs = max(n_yr_obs_r)) %>% 
  group_by(strat,stratum) %>% 
  summarise(n_routes = n(),
            mean_y_rt_obs = mean(mean_y_rt_obs),
            max_y_rt_obs = max(mean_y_rt_obs),
            mean_y_rt = mean(n_rt_obs),
            max_y_rt = max(n_rt_obs),
            .groups = "drop") %>% 
  distinct()

tr_f <- function(x){
  t <- (exp(x)-1)*100
}



trendst <- posterior_samples(fit = slope_stanfit,
                             parm = "beta",
                             dims = "strat") %>%
  posterior_sums(.,
                 dims = "strat")%>% 
  left_join(.,strat_df,by = "strat") %>% 
  mutate(trend = tr_f(mean),
         trend_lci = tr_f(lci),
         trend_uci = tr_f(uci),
         trend_se = tr_f(sd),
         trend_Wci = trend_uci-trend_lci)%>% 
  select(stratum,trend,trend_lci,trend_uci,trend_se,trend_Wci,n_routes)


trendst <- trendst %>% 
  mutate(version = spp)

save(list = "trendst",
     file = "data/Stan_slope_stratum_trends.RData")

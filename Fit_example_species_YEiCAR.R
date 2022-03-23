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

for(spp in c("YE_iCAR")){
out_base <- paste0(species_f,spp,firstYear,"_",lastYear)


# SPECIES LOOP ------------------------------------------------------------

output_dir <- "output"


sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")

    load(sp_file) 
output_dir <- "output"
out_base <- paste0(species_f,spp,firstYear,"_",lastYear)


  mod.file = "models/slopeYE_iCAR_route.stan"
    
    init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts,0,0.1),
                                 alpha_raw = rnorm(stan_data$nroutes,0,0.1),
                                 ALPHA = 0,
                                 BETA = 0,
                                 eta = 0,
                                 obs_raw = rnorm(stan_data$nobservers,0,0.1),
                                 sdnoise = 0.2,
                                 sd_year = 0.1,
                                 year_effect_raw = matrix(data = rnorm(n = stan_data$nroutes * stan_data$nyears,
                                                                       0,0.01),
                                                          ncol = stan_data$nyears,
                                                          nrow = stan_data$nroutes),
                                 sdobs = 0.1,
                                 sdbeta_space = runif(1,0.01,0.1),
                                 beta_raw_space = rnorm(stan_data$nroutes,0,0.01))} 
    
    
  



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


}



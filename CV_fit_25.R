## 1-step ahead, cross-validation of three route-level trend models for the BBS
setwd("C:/GitHub/BBS_iCAR_route_trends/")
library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(sf)
library(spdep)
library(ggforce)
source("functions/neighbours_define.R") ## function to define neighbourhood relationships
source("functions/prepare-jags-data-alt.R") ## small alteration of the bbsBayes function
source("functions/get_basemap_function.R") ## loads one of the bbsBayes strata maps
source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output
source("functions/initial_value_functions.R")
## changes captured in a commit on Nov 20, 2020


# load and stratify CASW data ---------------------------------------------
#species = "Pacific Wren"
#species = "Barred Owl"
strat = "bbs_usgs"
model = "slope"
scope = "RangeWide"

firstYear = 2004
lastYear = 2019 # final year to consider

# select a minimum year for prediction (i.e., a route has to have data between 2004 and 2011 to be included)
# similar to "L" in Burkner et al 2020 (https://doi.org/10.1080/00949655.2020.1783262)
minimumYear = 2011 

load("Data/sp_sel.RData")


for(species in sp_sel){

  species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)
  species_f <- gsub(species_f,pattern = "'",replacement = "",fixed = T)
  

# CROSS-VALIDATION loop through the annual re-fitting --------------------------------------


for(sppn in c("iCAR","BYM","Non_spatial")){
  
  load(paste0("data/",species_f,"CV_base_data.RData"))
  
  output_dir <- "output"
  spp <- paste0("_",sppn,"_")
  
  

  predictions_save <- NULL
  

for(ynext in (minimumYear+1):lastYear){
  
  out_base <- paste0(species_f,spp,firstYear,"_",ynext,"_CV")
  
  sp_file <- paste0(output_dir,"/",out_base,".RData")
  
  # setting up the fitting data ------------------------------------------
  
  
  obs_df_fit <- full_obs_df[which(full_obs_df$r_year <= ynext-1),]
  
  stan_data <- list(count = obs_df_fit$count,
                    year = obs_df_fit$year,
                    route = obs_df_fit$routeF,
                    firstyr = obs_df_fit$firstyr,
                    observer = obs_df_fit$observer,
                    nobservers = max(obs_df_fit$observer),
                    nyears = max(obs_df_fit$year),
                    nroutes = nroutes,
                    ncounts = length(obs_df_fit$count),
                    fixedyear = floor(max(obs_df_fit$year)/2))

  if(spp != "_Non_spatial_"){
    
  stan_data[["N_edges"]] = car_stan_dat$N_edges
  stan_data[["node1"]] = car_stan_dat$node1
  stan_data[["node2"]] = car_stan_dat$node2
  }  
  
  # setting up the prediction data ------------------------------------------
  
  
  obs_df_predict <- full_obs_df[which(full_obs_df$r_year == ynext),]
  
  stan_data[["route_pred"]] <- obs_df_predict$routeF
  stan_data[["count_pred"]] <- obs_df_predict$count
  stan_data[["firstyr_pred"]] <- obs_df_predict$firstyr
  stan_data[["observer_pred"]] <- obs_df_predict$observer
  stan_data[["ncounts_pred"]] <- length(obs_df_predict$count)
  
  
  mod.file = paste0("models/slope",spp,"route_LFO_CV.stan")
  
  ## compile model
  slope_model <- cmdstan_model(mod.file)
  
  
  init_def <- init_func_list[[sppn]]
  
  
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
  


  

  
  log_lik_samples_full <- posterior_samples(fit = slope_stanfit,
                                            parm = "log_lik",
                                            dims = "i") 
  
  log_lik_samples <- log_lik_samples_full %>% 
    posterior_sums(.,quantiles = NULL,dims = "i") 
  names(log_lik_samples) <- paste0("log_lik_",names(log_lik_samples))
  
  
  E_pred_samples_full <- posterior_samples(fit = slope_stanfit,
                                           parm = "E_pred",
                                           dims = "i") 
  
  E_pred_samples <- E_pred_samples_full %>% 
    posterior_sums(.,quantiles = NULL,dims = "i") 
  names(E_pred_samples) <- paste0("E_pred_",names(E_pred_samples))
  
  
  obs_df_predict_out <- bind_cols(obs_df_predict,log_lik_samples)
  obs_df_predict_out <- bind_cols(obs_df_predict_out,E_pred_samples)
  obs_df_predict_out$species <- species
  obs_df_predict_out$model <- sppn
  obs_df_predict_out$base <- out_base
  
  
  
  predictions_save <- bind_rows(predictions_save,obs_df_predict_out)
  

  print(paste("Finished",sppn,ynext))
  
  
  save(list = c("predictions_save"),file = paste0("output/",species_f,spp,"_pred_save.RData"))
  
  
}
  
  
 




}
  
}#end species loop

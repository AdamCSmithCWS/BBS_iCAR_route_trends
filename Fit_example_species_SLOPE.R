## building a BYM route-level trend model for the BBS
setwd("C:/GitHub/BBS_iCAR_route_trends")
library(tidyverse)
library(bbsBayes)
# library(cmdstanr)
# # library(rstan)
# # rstan_options(auto_write = TRUE, javascript = FALSE)
# # library(shinystan)
# library(sf)
# library(spdep)
# # library(doParallel)
# # library(foreach)
#  library(ggforce)
# #library(tidybayes)
# #source("functions/mungeCARdata4stan.R")
# source("functions/neighbours_define.R") ## function to define neighbourhood relationships
# source("functions/prepare-jags-data-alt.R") ## small alteration of the bbsBayes function
# source("functions/get_basemap_function.R") ## loads one of the bbsBayes strata maps
# source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output
# ## changes captured in a commit on Nov 20, 2020


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


  spp <- "SLOPE"
  
#   out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
# 
# 
# # SPECIES LOOP ------------------------------------------------------------
# 
# output_dir <- "output"
# 
# 
# sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")
# 
#     load(sp_file) 
# output_dir <- "output"
# out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
# 
# stan_data$ncounts
# jags_data$ncounts

strat_data = stratify(by = strat)


jags_data = bbsBayes::prepare_data(strat_data = strat_data,
                                   species_to_run = species,
                                   model = model,
                                   #n_knots = 10,
                                   min_year = firstYear,
                                   max_year = lastYear,
                                   min_n_routes = 1)# spatial neighbourhood define --------------------------------------------


jags_mod = run_model(jags_data = jags_data,
                     parameters_to_save = c("nslope"),
                     model_file_path = "temp.txt",
                     parallel = TRUE,
                     n_chains = 3)



inds_n <- generate_indices(jags_mod = jags_mod,
                           jags_data = jags_data)


inds_slope <- generate_indices(jags_mod = jags_mod,
                           jags_data = jags_data,
                           alternate_n = "nslope")

t_slope <- generate_trends(inds_slope)
n_slope <- generate_trends(inds_n)

write.csv(t_slope,file = "data/bbsBayes_slope_trends_example.csv",
          row.names = FALSE)

write.csv(inds_slope$data_summary,file = "data/bbsBayes_slope_indices_example.csv",
          row.names = FALSE)




# model with no year-effects ----------------------------------------------




jags_mod = run_model(jags_data = jags_data,
                     parameters_to_save = c("nslope"),
                     model_file_path = "temp_slope.txt",
                     parallel = TRUE,
                     n_chains = 3)





inds_slope <- generate_indices(jags_mod = jags_mod,
                               jags_data = jags_data,
                               alternate_n = "nslope")

t_slope <- generate_trends(inds_slope)

write.csv(t_slope,file = "data/bbsBayes_slope_only_trends_example.csv",
          row.names = FALSE)

write.csv(inds_slope$data_summary,file = "data/bbsBayes_slope_only_indices_example.csv",
          row.names = FALSE)





# spatial patterns in mean trends among species ---------------------------


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

strat_data = stratify(by = strat)

firstYear = 2004
lastYear = 2019



allspecies.eng = strat_data$species_strat$english

species_list = allspecies.eng#[-which(allspecies.eng %in% species_list)]



# optional removing the non-Canadian data ------------------------------------------
Canada_only <- FALSE

if(Canada_only){
  scope = "Canada"
  names_strata <- get_composite_regions(strata_type = strat) ## bbsBayes function that gets a list of the strata names and their composite regions (provinces, BCRs, etc.)
  
  us_strata_remove <- names_strata[which(names_strata$national == "US"),"region"] # character vector of the strata names to remove from data
}else{
  scope = "RangeWide"
  us_strata_remove <- NULL
}




trends_out2 <- read.csv(file = paste0("output/combined_",firstYear,"_",lastYear,"_",scope,"_trends_and_intercepts2.csv"))


trends_out_space2 <- read.csv(file = paste0("output/combined_",firstYear,"_",lastYear,"_",scope,"_spatial_trends_and_intercepts2.csv"))





groups = c("Grassland.birds",
           "Forest.birds",
           "Other.birds",
           "Aerial.insectivores",
           "suburban",
           "other.wetland.birds")



# group loop --------------------------------------------------------------

g_sel = groups[1]


# species selection of trend data -----------------------------------------


socb_c = which(grepl(names(socb),pattern = g_sel))

sp_sel1 = socb[which(socb[,socb_c] == "Included in group"),"species"]

sp_sel2 = sp_sel1[which(sp_sel1 %in% sps_all)]


# covariate data ----------------------------------------------------------

load("data/fp_change.RData")

# non-breeding range change in footprint by species
sp_w_fp_ch_dat <- fp_ch_non_breeding$english
#some species don't have range data on footprint
sp_sel = sp_sel2[which(sp_sel2 %in% sp_w_fp_ch_dat)]


dat = all %>% filter(species %in% sp_sel)









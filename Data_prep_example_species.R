## building a BYM route-level trend model for the BBS

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

scope = "RangeWide"


species = "Blue-headed Vireo"
#species = "Dickcissel"

species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)

 
  sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")
 
 jags_data = prepare_jags_data(strat_data = strat_data,
                             species_to_run = species,
                             model = model,
                             #n_knots = 10,
                             min_year = firstYear,
                             max_year = lastYear,
                             min_n_routes = 1)# spatial neighbourhood define --------------------------------------------

 # strata map of one of the bbsBayes base maps
 # helps group and set boundaries for the route-level neighbours
  strata_map  <- get_basemap(strata_type = strat,
                         transform_laea = TRUE,
                         append_area_weights = FALSE)
 
 
realized_strata_map = filter(strata_map,ST_12 %in% unique(jags_data$strat_name))

# Spatial boundaries set up --------------------

# the iCAR (intrinsic Conditional AutoRegressive) spatial model uses neighbourhood structure
# to share information on abundance and trend (intercept and slope) among BBS routes
# 

strata_list <- data.frame(ST_12 = unique(jags_data$strat_name),
                          strat = unique(jags_data$strat))


realized_strata_map <- inner_join(realized_strata_map,strata_list, by = "ST_12")


strata_bounds <- st_union(realized_strata_map) #union to provide a simple border of the realised strata
strata_bounds_buf = st_buffer(strata_bounds,dist = 300000) #buffering the realised strata by 300km



jags_data[["routeF"]] <- as.integer(factor((jags_data$route)))

route_map = unique(data.frame(route = jags_data$route,
                              routeF = jags_data$routeF,
                              strat = jags_data$strat_name,
                              Latitude = jags_data$Latitude,
                              Longitude = jags_data$Longitude))


# reconcile duplicate spatial locations -----------------------------------
# adhoc way of separating different routes with the same starting coordinates
# this shifts the starting coordinates of teh duplicates by ~1.5km to the North East 
# ensures that the duplicates have a unique spatial location, but remain very close to
# their original location and retain the correct neighbourhood relationships
# these duplicates happen when a "new" route is established because some large proportion
# of the end of a route is changed, but the start-point remains the same
dups = which(duplicated(route_map[,c("Latitude","Longitude")]))
while(length(dups) > 0){
  route_map[dups,"Latitude"] <- route_map[dups,"Latitude"]+0.01 #=0.01 decimal degrees ~ 1km
  route_map[dups,"Longitude"] <- route_map[dups,"Longitude"]+0.01 #=0.01 decimal degrees ~ 1km
  dups = which(duplicated(route_map[,c("Latitude","Longitude")]))
  
}
dups = which(duplicated(route_map[,c("Latitude","Longitude")])) 
if(length(dups) > 0){stop(paste(spec,"ERROR - At least one duplicate route remains"))}


route_map = st_as_sf(route_map,coords = c("Longitude","Latitude"))
st_crs(route_map) <- 4269 #NAD83 commonly used by US federal agencies
#load strata map




route_map = st_transform(route_map,crs = st_crs(realized_strata_map))


## returns the adjacency data necessary for the stan model
## also exports maps and saved data objects to plot_dir
car_stan_dat <- neighbours_define(real_strata_map = route_map,
                  strat_link_fill = 100000,
                  plot_neighbours = TRUE,
                  species = species,
                  plot_dir = "data/",
                  plot_file = paste0("_",scope,"_route_maps"),
                  save_plot_data = TRUE,
                  voronoi = TRUE,
                  alt_strat = "routeF",
                  add_map = realized_strata_map)


 



stan_data = jags_data[c("ncounts",
                        #"nstrata",
                        #"nobservers",
                        "count",
                        #"strat",
                        #"obser",
                        "year",
                        "firstyr",
                        "fixedyear")]
stan_data[["nyears"]] <- max(jags_data$year)
stan_data[["observer"]] <- as.integer(factor((jags_data$ObsN)))
stan_data[["nobservers"]] <- max(stan_data$observer)



stan_data[["N_edges"]] = car_stan_dat$N_edges
stan_data[["node1"]] = car_stan_dat$node1
stan_data[["node2"]] = car_stan_dat$node2
stan_data[["route"]] = jags_data$routeF
stan_data[["nroutes"]] = max(jags_data$routeF)


if(car_stan_dat$N != stan_data[["nroutes"]]){stop("Some routes are missing from adjacency matrix")}

save(list = c("stan_data",
              "jags_data",
              "route_map",
              "realized_strata_map",
              "firstYear"),
     file = sp_file)




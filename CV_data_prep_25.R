## 1-step ahead, cross-validation of three route-level trend models for the BBS

library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(sf)
library(spdep)
library(ggforce)
source("functions/neighbours_define_alt.R") ## function to define neighbourhood relationships
source("functions/prepare-jags-data-alt.R") ## small alteration of the bbsBayes function
source("functions/get_basemap_function.R") ## loads one of the bbsBayes strata maps
source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output
## changes captured in a commit on Nov 20, 2020


# load and stratify CASW data ---------------------------------------------
#species = "Pacific Wren"
#species = "Barred Owl"
strat = "bbs_usgs"
model = "slope"
scope = "RangeWide"
strat_data = stratify(by = strat)

firstYear = 2004
lastYear = 2019 # final year to consider

# select a minimum year for prediction (i.e., a route has to have data between 2004 and 2011 to be included)
# similar to "L" in Burkner et al 2020 (https://doi.org/10.1080/00949655.2020.1783262)
minimumYear = 2011 


species_list <- strat_data$species_strat$english

#series of if statements that skip analysing hybrids, composite species groups, etc.
if(grepl(pattern = "hybrid",x = species)){next}
if(grepl(pattern = "unid.",x = species)){next}
if(grepl(pattern = " x ",x = species)){next}
if(grepl(pattern = "/",fixed = TRUE,x = species)){next}
if(substr(x = species,1,1) == "("){next}


pred_save_allsp <- NULL

sp_w_trends <- list.files(path = "Figures/images/")
sp_w_trends <- gsub(gsub(sp_w_trends,pattern = "_Trends_2004.png",replacement = ""),pattern = "_",replacement = " ")

rselsp <- sample(sp_w_trends,50)

selSpecies <- unique(c(selSpecies,rselsp))


for(species in selSpecies[c(54:40)]){
  
  species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)
  species_f <- gsub(species_f,pattern = "'",replacement = "",fixed = T)
  
  # if(file.exists(sp_file)){next}
  #   

  


jags_data_inc = prepare_jags_data(strat_data = strat_data,
                                      species_to_run = species,
                                      model = model,
                                      #n_knots = 10,
                                      min_year = firstYear,
                                      max_year = minimumYear,
                                      min_n_routes = 1) # this final argument removes all data from the US


# strata map of one of the bbsBayes base maps
# helps group and set boundaries for the route-level neighbours
strata_map  <- get_basemap(strata_type = strat,
                           transform_laea = TRUE,
                           append_area_weights = FALSE)


realized_strata_map = filter(strata_map,ST_12 %in% unique(jags_data_inc$strat_name))

# Spatial boundaries set up --------------------

# the iCAR (intrinsic Conditional AutoRegressive) spatial model uses neighbourhood structure
# to share information on abundance and trend (intercept and slope) among BBS routes
# 

strata_list <- data.frame(ST_12 = unique(jags_data_inc$strat_name),
                          strat = unique(jags_data_inc$strat))


realized_strata_map <- inner_join(realized_strata_map,strata_list, by = "ST_12")


strata_bounds <- st_union(realized_strata_map) #union to provide a simple border of the realised strata
strata_bounds_buf = st_buffer(strata_bounds,dist = 300000) #buffering the realised strata by 300km



jags_data_inc[["routeF"]] <- as.integer(factor((jags_data_inc$route)))

route_map = unique(data.frame(route = jags_data_inc$route,
                              routeF = jags_data_inc$routeF,
                              strat = jags_data_inc$strat_name,
                              Latitude = jags_data_inc$Latitude,
                              Longitude = jags_data_inc$Longitude))


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
                                  plot_dir = "Figures/",
                                  plot_file = paste0("_CV_",scope,"_route_maps"),
                                  save_plot_data = TRUE,
                                  voronoi = TRUE,
                                  alt_strat = "routeF",
                                  add_map = realized_strata_map)






# Remove the remaining routes from the bbsBayes stratified data -----------
routes_inc <- route_map %>% 
  data.frame() %>% 
  select(route,routeF,strat) 


strat_data_reduced <- strat_data
strat_data_reduced$route_strat <- strat_data_reduced$route_strat[which(strat_data_reduced$route_strat$rt.uni %in% routes_inc$route),]
strat_data_reduced$bird_strat <- strat_data_reduced$bird_strat[which(strat_data_reduced$bird_strat$rt.uni %in% routes_inc$route),]


jags_data_red_allyears <- prepare_jags_data(strat_data = strat_data_reduced,
                                            species_to_run = species,
                                            model = model,
                                            #n_knots = 10,
                                            min_year = firstYear,
                                            max_year = lastYear,
                                            min_n_routes = 1)

routes_inc_red <- unique(data.frame(strat_name = jags_data_red_allyears$strat_name,
                                    route = jags_data_red_allyears$route))



if(!any(routes_inc_red$route %in% routes_inc$route)){
  stop(paste("Some routes in recent data are missing from original data"))
}

full_obs_df <- data.frame(count = jags_data_red_allyears$count,
                          r_year = jags_data_red_allyears$r_year,
                          year = jags_data_red_allyears$year,
                          firstyr = jags_data_red_allyears$firstyr,
                          route = jags_data_red_allyears$route,
                          ObsN = jags_data_red_allyears$ObsN) %>% 
  left_join(.,routes_inc, by = c("route")) %>% 
  arrange(year)

#identify the order in which the observers show up in the dataset
obs_sort <- unique(full_obs_df$ObsN)
# create an observer index that is constant across time
full_obs_df$observer <- as.integer(factor(full_obs_df$ObsN,levels = obs_sort,ordered = TRUE))

nroutes <- max(routes_inc$routeF)







save(list = c("full_obs_df","route_map","routes_inc","car_stan_dat",
              "nroutes"),
     file = paste0("data/",species_f,"CV_base_data.RData"))




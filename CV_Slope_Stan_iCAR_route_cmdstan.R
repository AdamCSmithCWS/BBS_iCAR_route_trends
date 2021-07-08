## 1-step ahead, cross-validation of a BYM route-level trend model for the BBS
## example using two species. 

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
lastYear = 2019 # final year to consider

# select a minimum year for prediction (i.e., a route has to have data between 2004 and 2011 to be included)
# similar to "L" in Burkner et al 2020 (https://doi.org/10.1080/00949655.2020.1783262)
minimumYear = 2011 


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


# SPECIES LOOP ------------------------------------------------------------

output_dir <- "G:/BBS_iCAR_route_trends/output"
#output_dir_simp <- "G:/BBS_iCAR_route_trends/output"

#species = species_list[748]

# for(species in rev(allspecies.eng)){
#   


species <- "Tennessee Warbler"



species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)

# if(file.exists(sp_file)){next}
#   
#series of if statements that skip analysing hybrids, composite species groups, etc.
if(grepl(pattern = "hybrid",x = species)){next}
if(grepl(pattern = " x ",x = species)){next}
if(grepl(pattern = "/",fixed = TRUE,x = species)){next}
if(substr(x = species,1,1) == "("){next}



## identify the routes that can be included
## routes that have data between firstYear and mimimumYear
jags_data_inc = try(prepare_jags_data(strat_data = strat_data,
                                      species_to_run = species,
                                      model = model,
                                      #n_knots = 10,
                                      min_year = firstYear,
                                      max_year = minimumYear,
                                      min_n_routes = 1,
                                      strata_rem = us_strata_remove),silent = TRUE) # this final argument removes all data from the US



### now just hte Canadian data remain.
if(class(jags_data_inc) == "try-error"){next}
if(jags_data_inc$ncounts < 500){next}

# spatial neighbourhood define --------------------------------------------

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
                                  plot_dir = "route_maps/",
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











# CROSS-VALIDATION loop through the annual re-fitting --------------------------------------

predictions_save_CAR <- NULL
predictions_save_NonCAR <- NULL


for(ynext in (minimumYear+1):lastYear){
  
  sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",ynext,"_CV_iCAR.RData")
  
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
  stan_data[["N_edges"]] = car_stan_dat$N_edges
  stan_data[["node1"]] = car_stan_dat$node1
  stan_data[["node2"]] = car_stan_dat$node2
  
  
  # setting up the prediction data ------------------------------------------
  
  
  obs_df_predict <- full_obs_df[which(full_obs_df$r_year == ynext),]
  
  stan_data[["route_pred"]] <- obs_df_predict$routeF
  stan_data[["count_pred"]] <- obs_df_predict$count
  stan_data[["firstyr_pred"]] <- obs_df_predict$firstyr
  stan_data[["observer_pred"]] <- obs_df_predict$observer
  stan_data[["ncounts_pred"]] <- length(obs_df_predict$count)
  
  
  mod.file = "models/slope_iCAR_route2_LFO_CV.stan"
  
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
                               sdbeta_space = runif(1,0.01,0.1),
                               sdbeta_rand = runif(1,0.01,0.1),
                               beta_raw_space = rnorm(stan_data$nroutes,0,0.01),
                               beta_raw_rand = rnorm(stan_data$nroutes,0,0.01))}
  
  
  slope_stanfit <- slope_model$sample(
    data=stan_data,
    refresh=100,
    chains=3, iter_sampling=1000,
    iter_warmup=1000,
    parallel_chains = 3,
    #pars = parms,
    adapt_delta = 0.8,
    max_treedepth = 14,
    seed = 123,
    init = init_def)
  
  out_base <- paste0(species_f,"_",scope,"_",firstYear,"_CV_CAR_",ynext)
  
  # export to csv and read in as rstan --------------------------------------
  slope_stanfit$save_output_files(dir = output_dir,
                                  basename = out_base,
                                  random = FALSE,
                                  timestamp = FALSE)
  
  csv_files <- dir(output_dir,pattern = out_base,full.names = TRUE)
  
  
  shiny_explore <- FALSE
  if(shiny_explore){
    sl_rstan <- rstan::read_stan_csv(csv_files)
    launch_shinystan(as.shinystan(sl_rstan))
    
    loo_stan = loo(sl_rstan)
  }
  
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
  
  predictions_save_CAR <- bind_rows(predictions_save_CAR,obs_df_predict_out)
  
  
  
  
  
  # nonSpatial model fit ----------------------------------------------------
  
  mod.file.non = "models/slope_noniCAR_route2_LFO_CV.stan"
  
  ## compile model
  slope_model_non <- cmdstan_model(mod.file.non)
  
  init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts,0,0.1),
                               alpha_raw = rnorm(stan_data$nroutes,0,0.1),
                               ALPHA = 0,
                               BETA = 0,
                               eta = 0,
                               obs_raw = rnorm(stan_data$nobservers,0,0.1),
                               sdnoise = 0.2,
                               sdobs = 0.1,
                               #sdbeta_space = runif(1,0.01,0.1),
                               sdbeta_rand = runif(1,0.01,0.1),
                               #beta_raw_space = rnorm(stan_data$nroutes,0,0.01),
                               beta_raw_rand = rnorm(stan_data$nroutes,0,0.01))}
  
  
  slope_stanfit_non <- slope_model_non$sample(
    data=stan_data,
    refresh=100,
    chains=3, iter_sampling=1000,
    iter_warmup=1000,
    parallel_chains = 3,
    #pars = parms,
    adapt_delta = 0.8,
    max_treedepth = 14,
    seed = 123,
    init = init_def)
  
  
  
  out_base_non <- paste0(species_f,"_",scope,"_",firstYear,"_CV_non_",ynext)
  
  # export to csv and read in as rstan --------------------------------------
  slope_stanfit_non$save_output_files(dir = output_dir,
                                      basename = out_base_non,
                                      random = FALSE,
                                      timestamp = FALSE)
  
  csv_files <- dir(output_dir,pattern = out_base_non,full.names = TRUE)
  
  
  
  log_lik_samples_full_non <- posterior_samples(fit = slope_stanfit_non,
                                                parm = "log_lik",
                                                dims = "i") 
  
  log_lik_samples <- log_lik_samples_full_non %>% 
    posterior_sums(.,quantiles = NULL,dims = "i") 
  names(log_lik_samples) <- paste0("log_lik_",names(log_lik_samples))
  
  
  E_pred_samples_full_non <- posterior_samples(fit = slope_stanfit_non,
                                               parm = "E_pred",
                                               dims = "i") 
  
  E_pred_samples <- log_lik_samples_full_non %>% 
    posterior_sums(.,quantiles = NULL,dims = "i") 
  names(E_pred_samples) <- paste0("E_pred_",names(E_pred_samples))
  
  
  obs_df_predict_out <- bind_cols(obs_df_predict,log_lik_samples)
  obs_df_predict_out <- bind_cols(obs_df_predict_out,E_pred_samples)
  
  predictions_save_NonCAR <- bind_rows(predictions_save_NonCAR,obs_df_predict_out)
  
  
  
  
  save(list = c("log_lik_samples_full",
                "E_pred_samples_full",
                "slope_stanfit",
                "log_lik_samples_full_non",
                "E_pred_samples_full_non",
                "slope_stanfit_non",
                "stan_data",
                "obs_df_predict",
                "predictions_save_NonCAR",
                "predictions_save_CAR",
                "ynext",
                "sp_file",
                "species_f",
                "output_dir",
                "out_base_non",
                "out_base"),
       file = sp_file)
  
  print(ynext)
}





#stopCluster(cl = cluster)








# post loop analysis ------------------------------------------------------


# 
# launch_shinystan(slope_stanfit) 
# 
# 
# library(loo)
# library(tidyverse)
# 
# log_lik_1 <- extract_log_lik(slope_stanfit, merge_chains = FALSE)
# r_eff <- relative_eff(exp(log_lik_1), cores = 10)
# loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 10)
# print(loo_1)
# 
# doy = ((jags_data$month-4)*30+jags_data$day)
# plot(loo_1$pointwise[,"influence_pareto_k"],log(stan_data$count+1))
# plot(loo_1$pointwise[,"influence_pareto_k"],doy)
# plot(doy,log(stan_data$count+1))
# 
# 
# 
# loo2 = data.frame(loo_1$pointwise)
# 
# loo2$flag = cut(loo2$influence_pareto_k,breaks = c(0,0.5,0.7,1,Inf))
# dts = data.frame(count = stan_data$count,
#                  obser = stan_data$obser,
#                  route = stan_data$route,
#                  year = stan_data$year)
# loo2 = cbind(loo2,dts)
# 
# plot(log(loo2$count+1),loo2$influence_pareto_k)
# 
# obserk = loo2 %>% group_by(obser) %>% 
#   summarise(n = log(n()),
#             mean_k = mean(influence_pareto_k),
#             max_k = max(influence_pareto_k),
#             sd_k = sd(influence_pareto_k),
#             mean_looic = mean(looic),
#             mean_ploo = mean(p_loo))
# plot(obserk$n,obserk$max_k)
# plot(obserk$n,obserk$mean_k)
# plot(obserk$n,obserk$sd_k)
# plot(obserk$n,obserk$mean_looic)
# plot(obserk$n,obserk$mean_ploo)
# 
# 
# yeark = loo2 %>% group_by(year) %>% 
#   summarise(n = n(),
#             mean_k = mean(influence_pareto_k),
#             q90 = quantile(influence_pareto_k,0.9),
#             max_k = max(influence_pareto_k),
#             sd_k = sd(influence_pareto_k),
#             route = mean(route),
#             sd = sd(route))
# plot(yeark$year,yeark$max_k)
# plot(yeark$year,yeark$mean_k)
# plot(yeark$year,yeark$sd_k)
# plot(yeark$year,yeark$q90)
# 
# routek = loo2 %>% group_by(route) %>% 
#   summarise(n = n(),
#             mean_k = mean(influence_pareto_k),
#             q90_k = quantile(influence_pareto_k,0.9),
#             max_k = max(influence_pareto_k),
#             sd_k = sd(influence_pareto_k),
#             route = mean(route),
#             sd = sd(route))
# plot(routek$route,routek$max_k)
# plot(routek$n,routek$mean_k)
# 
# plot(routek$route,routek$mean_k)
# plot(routek$route,routek$sd_k)
# plot(routek$route,routek$q90_k)
# 
# 


# PLOTTING and trend output -----------------------------------------------

# library(tidybayes)


route_trajectories <- FALSE #set to FALSE to speed up mapping

maps = vector(mode = "list",length = 400)
maps2 = vector(mode = "list",length = 400)

maps3 = vector(mode = "list",length = 400)

maps_rand = vector(mode = "list",length = 400)
maps_space = vector(mode = "list",length = 400)

trends_out <- NULL
trends_out_space <- NULL
trends_out_rand <- NULL
sdbeta_dif <- NULL
sdbeta_space_rand <- NULL

jj <- 0

LC = 0.05
UC = 0.95

output_dir <- "G:/BBS_iCAR_route_trends/output"

strata_map  <- get_basemap(strata_type = strat,
                           transform_laea = TRUE,
                           append_area_weights = FALSE)




for(species in rev(allspecies.eng)){
  
  
  species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)
  
  sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")
  
  
  #sp_file <- paste0("output/",species,"Canadian_",firstYear,"_",lastYear,"_slope_route_iCAR2.RData")
  if(file.exists(sp_file)){
    
    load(sp_file)
    # if(species == "Northern Cardinal"){next
    #   
    # #csv_files <- dir(output_dir,pattern = out_base,full.names = TRUE)
    # }
    if(length(csv_files) == 0 | length(csv_files > 3)){
      csv_files = paste0(output_dir,"/",
                         species_f,"_",scope,"_",firstYear,
                         "-",1:3,".csv")
    }
    ### may be removed after re-running     launch_shinystan(slope_stanfit)
    sl_rstan <- rstan::read_stan_csv(csv_files)
    #launch_shinystan(as.shinystan(sl_rstan))
    
    #loo_stan = loo(sl_rstan)
    
    jj <- jj+1
    # laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system
    # 
    # locat = system.file("maps",
    #                     package = "bbsBayes")
    # map.file = "BBS_USGS_strata"
    # 
    # strata_map = read_sf(dsn = locat,
    #                      layer = map.file)
    # strata_map = st_transform(strata_map,crs = laea)
    # 
    # realized_strata_map = filter(strata_map,ST_12 %in% unique(jags_data$strat_name))
    # 
    # strata_list <- data.frame(ST_12 = unique(jags_data$strat_name),
    #                           strat = unique(jags_data$strat))
    # 
    # 
    # realized_strata_map <- inner_join(realized_strata_map,strata_list, by = "ST_12")
    # 
    
    ####
    # add trend and abundance ----------------------------------------
    
    beta_samples = posterior_samples(sl_rstan,"beta",
                                     dims = "s")
    # beta_samples2 = posterior_samples(slope_stanfit,"beta",
    #                                  dims = "s")
    
    slopes = beta_samples %>% group_by(s) %>% 
      summarise(b = mean(.value),
                lci = quantile(.value,LC),
                uci = quantile(.value,UC),
                sd = sd(.value),
                prec = 1/var(.value),
                trend = mean((exp(.value)-1)*100),
                lci_trend = quantile((exp(.value)-1)*100,LC),
                uci_trend = quantile((exp(.value)-1)*100,UC),
                .groups = "keep")
    
    alpha_samples = posterior_samples(sl_rstan,"alpha",
                                      dims = "s")
    interc = alpha_samples %>% group_by(s) %>% 
      summarise(abund = mean(exp(.value)),
                lci_i = quantile(exp(.value),LC),
                uci_i = quantile(exp(.value),UC),
                sd_i = sd(exp(.value)),
                prec_i = 1/var(.value),
                .groups = "keep")
    
    #plot(log(interc$i),slopes$b)
    slops_int = inner_join(slopes,interc,by = "s")
    slops_int$routeF = slops_int$s
    
    
    
    # random effect plus mean component of slope ----------------------------------------
    
    # BETA_samples = posterior_samples(sl_rstan,BETA) %>% 
    #   rename(BETA = .value) %>% 
    #   ungroup() %>% 
    #   select(BETA,.draw)
    # 
    # beta_rand_samples = posterior_samples(sl_rstan,beta_rand[s]) %>% 
    #   rename(beta_rand = .value) %>% 
    #   ungroup() %>% 
    #   select(beta_rand,.draw,s)
    # 
    # beta_rand_samples <- inner_join(beta_rand_samples,BETA_samples,by = c(".draw"))
    # 
    # slopes_rand_full = beta_rand_samples %>% group_by(s) %>% 
    #   summarise(b = mean(beta_rand + BETA),
    #             lci = quantile(beta_rand + BETA,LC),
    #             uci = quantile(beta_rand + BETA,UC),
    #             sd = sd(beta_rand + BETA),
    #             prec = 1/var(beta_rand + BETA),
    #             .groups = "keep")
    # 
    # slopes_rand_full_int = inner_join(slopes_rand_full,interc,by = "s")
    # slopes_rand_full_int$routeF = slopes_rand_full_int$s
    
    
    beta_rand_samples = posterior_samples(sl_rstan,
                                          "beta_rand",
                                          dims = "s")
    
    slopes_rand = beta_rand_samples %>% group_by(s) %>% 
      summarise(b = mean(.value),
                lci = quantile(.value,LC),
                uci = quantile(.value,UC),
                sd = sd(.value),
                prec = 1/var(.value),
                trend = mean((exp(.value)-1)*100),
                lci_trend = quantile((exp(.value)-1)*100,LC),
                uci_trend = quantile((exp(.value)-1)*100,UC),
                .groups = "keep")
    
    slops_rand_int = inner_join(slopes_rand,interc,by = "s")
    slops_rand_int$routeF = slops_rand_int$s
    
    
    # spatial component of slope ----------------------------------------
    
    
    beta_space_samples = posterior_samples(sl_rstan,"beta_space",
                                           dims = "s")
    
    slopes_space = beta_space_samples %>% group_by(s) %>% 
      summarise(b = mean(.value),
                lci = quantile(.value,LC),
                uci = quantile(.value,UC),
                sd = sd(.value),
                prec = 1/var(.value),
                trend = mean((exp(.value)-1)*100),
                lci_trend = quantile((exp(.value)-1)*100,LC),
                uci_trend = quantile((exp(.value)-1)*100,UC),
                .groups = "keep")
    
    slops_space_int = inner_join(slopes_space,interc,by = "s")
    slops_space_int$routeF = slops_space_int$s
    
    
    # Compare spatial and random variation ------------------------------------
    sdbeta_rand_tmp_samples <- posterior_samples(sl_rstan,
                                                 "sdbeta_rand")
    sdbeta_space_tmp_samples <- posterior_samples(sl_rstan,
                                                  "sdbeta_space")
    
    sdbeta_space_rand_tmp_samples <- bind_rows(sdbeta_rand_tmp_samples,
                                               sdbeta_space_tmp_samples)
    
    
    sdbeta_space_rand_tmp <- sdbeta_space_rand_tmp_samples %>% 
      group_by(.variable) %>%
      summarise(mean = mean((.value)),
                lci = quantile((.value),LC),
                uci = quantile((.value),UC),
                sd = sd((.value)),
                .groups = "keep") %>% 
      mutate(species = species)
    #combines all species estimates
    sdbeta_space_rand <- bind_rows(sdbeta_space_rand,sdbeta_space_rand_tmp)
    
    
    
    # difference rand-spatial -------------------------------------------------
    
    sdbeta_space_tmp_samples <- sdbeta_space_tmp_samples %>% 
      rename(sd_space = .value) %>% 
      ungroup() %>% 
      select(-.variable)
    
    sdbeta_rand_tmp_samples <- sdbeta_rand_tmp_samples %>% 
      rename(sd_rand = .value)%>% 
      ungroup() %>% 
      select(-.variable)
    
    sdbeta_tmp_samples <- inner_join(sdbeta_rand_tmp_samples,sdbeta_space_tmp_samples)
    
    sdbeta_tmp_dif <- sdbeta_tmp_samples %>% 
      group_by(.draw) %>%
      summarise(dif = sd_rand-sd_space) %>% 
      ungroup() %>% 
      summarise(mean = mean((dif)),
                lci = quantile((dif),LC),
                uci = quantile((dif),UC),
                sd = sd((dif))) %>% 
      mutate(species = species)
    
    sdbeta_dif <- bind_rows(sdbeta_dif,sdbeta_tmp_dif)
    
    
    
    # Route-level trajectories ------------------------------------------------
    if(route_trajectories){
      sdnoise_samples = posterior_samples(sl_rstan,"sdnoise")%>% 
        ungroup() %>% 
        select(.draw,.value) %>% 
        rename(sdnoise = .value)
      
      sdobs_samples = posterior_samples(sl_rstan,"sdobs")%>% 
        ungroup() %>% 
        select(.draw,.value) %>% 
        rename(sdobs = .value)
      
      
      beta_samples <- beta_samples %>% 
        ungroup() %>% 
        select(s,.draw,.value) %>% 
        rename(beta = .value)
      
      alpha_samples <- alpha_samples %>% 
        ungroup() %>% 
        select(s,.draw,.value) %>% 
        rename(alpha = .value)
      
      ab_samples <- inner_join(beta_samples,alpha_samples)
      
      ab_samples <- ab_samples %>% 
        left_join(.,sdnoise_samples,by = ".draw") %>% 
        left_join(.,sdobs_samples,by = ".draw")
      
      
      nyears = stan_data$nyears
      fixedyear = stan_data$fixedyear
      YEARS = c(min(jags_data$r_year):max(jags_data$r_year))
      
      if(length(YEARS) != nyears){stop("years don't match YEARS =",length(YEARS),"nyears =",nyears)}
      
      ind_fxn = function(a,b,sdn,sdob,y,fy){
        i = exp(a + b*(y-fy) + (0.5*(sdn^2))+ (0.5*(sdob^2)))
        return(i)
      }
      
      ### this could be simplified to just estimate the start and end-years
      i_samples = NULL
      for(yr in 1:nyears){
        i_t = ab_samples %>% 
          mutate(i = ind_fxn(alpha,beta,sdnoise,sdobs,yr,fixedyear),
                 y = yr,
                 year = YEARS[yr])
        i_samples <- bind_rows(i_samples,i_t)
      }
      
      ### this could be tweaked to sum across all routes in the original strata
      ### just join to the strata-route dataframe - route_map
      ### then add the non-zero-weights for the strata
      ### then add the area-weights for the strata
      ### and change the group-by value
      indices = i_samples %>% group_by(s,y,year) %>% 
        summarise(index = mean(i),
                  lci_i = quantile(i,LC),
                  uci_i = quantile(i,UC),
                  sd_i = sd(i),
                  .groups = "keep")
      
      raw = data.frame(s = stan_data$route,
                       y = stan_data$year,
                       count = stan_data$count,
                       obs = stan_data$observer)
      indices = left_join(indices,raw,by = c("y","s"))
      
      rts = route_map %>% tibble() %>% 
        select(route,routeF,strat) %>% 
        mutate(s = routeF) 
      
      
      indices = left_join(indices,rts,by = "s")
      indices$obs <- factor(indices$obs)
      nroutes = stan_data$nroutes
      
      # setting up the plot dimensions
      npg = ceiling(nroutes/9)
      ncl = 3
      nrw = 3
      if(npg*9-nroutes < 3){
        nrw = 2
        npg = ceiling(nroutes/6) 
        if(npg*6-nroutes < 3){
          ncl = 2
          npg = ceiling(nroutes/4)  
        }
      }
      #### 
      pdf(paste0("trajectories/",species,"_route_trajectories2.pdf"),
          width = 11,
          height = 8.5)
      
      for(j in 1:npg){
        traj = ggplot(data = indices,aes(x = year,y = index,colour = strat))+
          geom_ribbon(aes(ymin = lci_i,ymax = uci_i),alpha = 0.4,fill = grey(0.5))+
          geom_line()+
          geom_point(aes(x = year,y = count, colour = obs), fill = grey(0.5),alpha = 0.5,inherit.aes = FALSE)+
          facet_wrap_paginate(~ strat+route,scales = "free",ncol = ncl,nrow = nrw,page = j)+
          theme(legend.position = "none")
        try(print(traj),silent = TRUE)
      }
      dev.off()
      
      
      
    }
    
    # connect trends to original route names ----------------------------------
    
    route_map_out = left_join(route_map,slops_int,by = "routeF")
    route_map_out$species <- species
    
    trends_out <- bind_rows(trends_out,route_map_out)
    
    
    
    route_map_out_rand = left_join(route_map,slops_rand_int,by = "routeF")
    route_map_out_rand$species <- species
    
    trends_out_rand <- bind_rows(trends_out_rand,route_map_out_rand)
    
    # slopes_rand_full_int
    # route_map_out_rand = left_join(route_map,slopes_rand_full_int,by = "routeF")
    # route_map_out_rand$species <- species
    # 
    # trends_out_rand <- bind_rows(trends_out_rand,route_map_out_rand)
    
    
    
    route_map_out_space = left_join(route_map,slops_space_int,by = "routeF")
    route_map_out_space$species <- species
    
    trends_out_space <- bind_rows(trends_out_space,route_map_out_space)
    
    
    
    ### setting up boundaries for plots
    # load(paste0("route_maps/",species_f,"_route_data.RData"))
    
    strata_bounds <- st_union(realized_strata_map) #union to provide a simple border of the realised strata
    bb = st_bbox(strata_bounds)
    xlms = as.numeric(c(bb$xmin,bb$xmax))
    ylms = as.numeric(c(bb$ymin,bb$ymax))
    
    
    
    
    # add mapping of trends ---------------------------------------------------
    
    plot_trend <- TRUE #set to false to plot the slopes
    if(plot_trend){
      breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
      lgnd_head <- "Trend\n"
      trend_title <- "trends"
      labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
      labls = paste0(labls, " %/year")
      route_map_out$Tplot <- cut(route_map_out$trend,breaks = c(-Inf, breaks, Inf),labels = labls)
      route_map_out <- route_map_out %>% 
        mutate(h_ci = (uci_trend-lci_trend)/2)
      
      route_map_out_space$Tplot <- cut(route_map_out_space$trend,breaks = c(-Inf, breaks, Inf),labels = labls)
      route_map_out_rand$Tplot <- cut(route_map_out_rand$trend,breaks = c(-Inf, breaks, Inf),labels = labls)
      
      
    }else{
      breaks <- c(-0.07, -0.04, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.04, 0.07)
      lgnd_head <- "slope\n"
      trend_title <- "trend-slopes"
      labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
      labls = paste0(labls, " slope")
      route_map_out$Tplot <- cut(route_map_out$b,breaks = c(-Inf, breaks, Inf),labels = labls)
      route_map_out <- route_map_out %>% 
        mutate(h_ci = (uci-lci)/2)
      
      route_map_out_space$Tplot <- cut(route_map_out_space$b,breaks = c(-Inf, breaks, Inf),labels = labls)
      route_map_out_rand$Tplot <- cut(route_map_out_rand$b,breaks = c(-Inf, breaks, Inf),labels = labls)
      
    }
    map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                     "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
    names(map_palette) <- labls
    
    
    tmap = ggplot(route_map_out)+
      #geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
      geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
      geom_sf(aes(colour = Tplot,size = abund))+
      scale_size_continuous(range = c(0.5,3),
                            name = "Mean abundance")+
      scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                          guide = guide_legend(reverse=TRUE),
                          name = paste0(lgnd_head,firstYear,"-",lastYear))+
      coord_sf(xlim = xlms,ylim = ylms)+
      labs(title = paste("DRAFT ",species,trend_title,"by BBS route"),
           subtitle = "Route-level trends from a spatial iCAR model, using Stan")
    
    maps[[jj]] <- tmap
    
    
    png(filename = paste0("Figures/images/",species_f,"_Trends_",firstYear,".png"),
        res = 600,
        width = 20,
        height = 15,
        units = "cm")
    print(tmap)
    dev.off()
    
    
    tmap2 = ggplot(route_map_out)+
      geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
      geom_sf(aes(colour = Tplot,size = abund))+
      scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                          guide = guide_legend(reverse=TRUE),
                          name = paste0(lgnd_head,firstYear,"-",lastYear))+
      coord_sf(xlim = xlms,ylim = ylms)+
      theme(legend.position = "none")+
      labs(title = paste(species))
    
    maps2[[jj]] <- tmap2
    
    
    
    tmap3 = ggplot(route_map_out)+
      geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
      geom_sf(aes(colour = Tplot,size = 1/h_ci))+
      scale_size_continuous(range = c(0.5,3))+
      scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                          guide = guide_legend(reverse=TRUE),
                          name = paste0(lgnd_head,firstYear,"-",lastYear))+
      coord_sf(xlim = xlms,ylim = ylms)+
      labs(title = paste(species))
    
    maps3[[jj]] <- tmap3
    
    
    
    
    tmap_space = ggplot(route_map_out_space)+
      geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
      geom_sf(aes(colour = Tplot,size = abund))+
      scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                          guide = guide_legend(reverse=TRUE),
                          name = paste0(lgnd_head,firstYear,"-",lastYear))+
      coord_sf(xlim = xlms,ylim = ylms)+
      theme(legend.position = "none")+
      labs(title = paste("spatial component"))
    maps_space[[jj]] <- tmap_space
    
    
    
    
    
    tmap_rand = ggplot(route_map_out_rand)+
      geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
      geom_sf(aes(colour = Tplot,size = abund))+
      scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                          guide = guide_legend(reverse=TRUE),
                          name = paste0(lgnd_head,firstYear,"-",lastYear))+
      coord_sf(xlim = xlms,ylim = ylms)+
      theme(legend.position = "none")+
      labs(title = paste("random component"))
    
    maps_rand[[jj]] <- tmap_rand
    
    
    print(species)
    # write.csv(route_map_out,
    #           file = paste0("output/",species," ",firstYear," ",lastYear,"_Canadian_trends_and_intercepts.csv"))
    
    
  }
}



# overall trend maps and trends -------------------------------------------


pdf(file = paste0("figures/Combined_",firstYear,"_",lastYear,"_",scope,"_trend_map_route2.pdf"),
    height = 8.5,
    width = 11)
for(j in 1:length(maps)){
  if(!is.null(maps[[j]])){print(maps[[j]])}
}
dev.off()



# comparison trend maps and trends -------------------------------------------

library(patchwork)

pdf(file = paste0("figures/Combined_",firstYear,"_",lastYear,"_",scope,"_trend_map_route2_by_half_CI.pdf"),
    height = 8.5,
    width = 11)
for(j in 1:length(maps3)){
  if(!is.null(maps3[[j]])){print(maps3[[j]])}
}
dev.off()





pdf(file = paste0("figures/Combined_",firstYear,"_",lastYear,"_",scope,"_space_all_trend_map_route2.pdf"),
    width = 8.5,
    height = 11)
for(j in 1:length(maps)){
  if(!is.null(maps[[j]])){
    print(maps2[[j]] /(maps_space[[j]]))
  }
}
dev.off()





# rename the trend output columns -----------------------------------------

trends_out2 = trends_out %>% 
  mutate(h_ci = (uci_trend-lci_trend)/2) %>% 
  select(.,
         species,route,strat,
         trend,lci_trend,uci_trend,h_ci,
         abund,lci_i,uci_i,
         b,lci,uci,sd) %>% 
  relocate(.,
           species,route,strat,
           trend,lci_trend,uci_trend,h_ci,
           abund,lci_i,uci_i,
           b,lci,uci,sd) %>% 
  rename(.,
         english_name = species,BBS_route = route, BBS_stratum = strat,
         Trend = trend,lci95_Trend = lci_trend,uci95_Trend = uci_trend,half_CI_width = h_ci,
         Mean_abundance = abund,lci95_Mean_abundance = lci_i,uci95_Mean_abundance = uci_i,
         slope = b,lci95_slope = lci,uci95_slope = uci,sd_slope = sd)


trends_out_space2 = trends_out_space %>% 
  mutate(h_ci = (uci_trend-lci_trend)/2) %>% 
  select(.,
         species,route,strat,
         trend,lci_trend,uci_trend,h_ci,
         abund,lci_i,uci_i,
         b,lci,uci,sd) %>% 
  relocate(.,
           species,route,strat,
           trend,lci_trend,uci_trend,h_ci,
           abund,lci_i,uci_i,
           b,lci,uci,sd) %>% 
  rename(.,
         english_name = species,BBS_route = route, BBS_stratum = strat,
         Trend = trend,lci95_Trend = lci_trend,uci95_Trend = uci_trend,half_CI_width = h_ci,
         Mean_abundance = abund,lci95_Mean_abundance = lci_i,uci95_Mean_abundance = uci_i,
         slope = b,lci95_slope = lci,uci95_slope = uci,sd_slope = sd)

trends_out_rand2 = trends_out_rand %>% 
  mutate(h_ci = (uci_trend-lci_trend)/2) %>% 
  select(.,
         species,route,strat,
         trend,lci_trend,uci_trend,h_ci,
         abund,lci_i,uci_i,
         b,lci,uci,sd) %>% 
  relocate(.,
           species,route,strat,
           trend,lci_trend,uci_trend,h_ci,
           abund,lci_i,uci_i,
           b,lci,uci,sd) %>% 
  rename(.,
         english_name = species,BBS_route = route, BBS_stratum = strat,
         Trend = trend,lci95_Trend = lci_trend,uci95_Trend = uci_trend,half_CI_width = h_ci,
         Mean_abundance = abund,lci95_Mean_abundance = lci_i,uci95_Mean_abundance = uci_i,
         slope = b,lci95_slope = lci,uci95_slope = uci,sd_slope = sd)


# Export the trend estimates ----------------------------------------------


write.csv(trends_out2,
          file = paste0("output/combined_",firstYear,"_",lastYear,"_",scope,"_trends_and_intercepts2.csv"))


write.csv(trends_out_space2,
          file = paste0("output/combined_",firstYear,"_",lastYear,"_",scope,"_spatial_trends_and_intercepts2.csv"))
write.csv(trends_out_rand2,
          file = paste0("output/combined_",firstYear,"_",lastYear,"_",scope,"_random_trends_and_intercepts2.csv"))



# graph the spatial and random variance comparison ------------------------

var_plot = ggplot()+
  geom_boxplot(data = sdbeta_space_rand,aes(x = .variable,y = mean))+
  theme_classic()

sdbeta_difs <- sdbeta_dif %>% 
  mutate(species = fct_reorder(species,mean))

var_dif_plot = ggplot(data = sdbeta_difs,aes(x = species,y = mean))+
  geom_point(aes(size = 1/(sd^2)))+
  geom_errorbar(aes(ymin = lci, ymax = uci),alpha = 0.3,width = 0)+
  geom_hline(yintercept = 0)+
  scale_size_continuous(range = c(0.5,2))+
  labs(title = "Difference in sd_beta (random - spatial)")+
  theme(axis.text = element_text(size = 5),
        legend.position = "none")+
  coord_flip()
pdf(file = paste0("figures/Combined_",firstYear,"_",lastYear,"_",scope,"_difference_beta_sd.pdf"),
    width = 8.5,
    height = 17)
print(var_dif_plot)
dev.off()

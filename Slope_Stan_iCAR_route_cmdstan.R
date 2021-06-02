## building a Stan version of the bbsBayes models

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


# SPECIES LOOP ------------------------------------------------------------

output_dir <- "G:/BBS_iCAR_route_trends/output"
#output_dir_simp <- "G:/BBS_iCAR_route_trends/output"

#species = species_list[748]

for(species in rev(allspecies.eng)){
  
  species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)
  
  sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")
  if(file.exists(sp_file)){next}
    
  #series of if statements that skip analysing hybrids, composite species groups, etc.
  if(grepl(pattern = "hybrid",x = species)){next}
  if(grepl(pattern = " x ",x = species)){next}
  if(grepl(pattern = "/",fixed = TRUE,x = species)){next}
  if(substr(x = species,1,1) == "("){next}
  
  

 jags_data = try(prepare_jags_data(strat_data = strat_data,
                             species_to_run = species,
                             model = model,
                             #n_knots = 10,
                             min_year = firstYear,
                             max_year = lastYear,
                             min_n_routes = 1,
                             strata_rem = us_strata_remove),silent = TRUE) # this final argument removes all data from the US
### now just hte Canadian data remain.
 if(class(jags_data) == "try-error"){next}
 if(jags_data$ncounts < 500){next}

# spatial neighbourhood define --------------------------------------------

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
                  plot_dir = "route_maps/",
                  plot_file = paste0("_",scope,"_route_maps.pdf"),
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

mod.file = "models/slope_iCAR_route2.stan"

# parms = c("sdnoise",
#           "sdobs",
#           "sdbeta_rand",
#           "sdbeta_space",
#           "sdalpha",
#           "BETA",
#           "ALPHA",
#           "beta",
#           "beta_rand",
#           "beta_space",
#           "alpha",
#           "eta",
#           "log_lik")

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
  refresh=25,
  chains=3, iter_sampling=1000,
  iter_warmup=1000,
  parallel_chains = 3,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 14,
  seed = 123,
  init = init_def)

out_base <- paste0(species_f,"_",scope,"_",firstYear)

# export to csv and read in as rstan --------------------------------------
slope_stanfit$save_output_files(dir = output_dir,
                                basename = out_base,
                                random = FALSE,
                                timestamp = FALSE)

csv_files <- dir(output_dir,pattern = out_base,full.names = TRUE)

#slope_stanfit$save_object(file = paste0(output_dir,"/",out_base,"_gamye_iCAR.RDS"))



shiny_explore <- FALSE
if(shiny_explore){
  sl_rstan <- rstan::read_stan_csv(csv_files)
  launch_shinystan(as.shinystan(sl_rstan))
  
  loo_stan = loo(sl_rstan)
}





# slope_model = stan_model(file=mod.file)
# 
# print(paste(firstYear,species))
# ## run sampler on model, data
# slope_stanfit <- sampling(slope_model,
#                                data=stan_data,
#                                verbose=TRUE, refresh=100,
#                                chains=4, iter=900,
#                                warmup=600,
#                                cores = 4,
#                                pars = parms,
#                                control = list(adapt_delta = 0.8,
#                                               max_treedepth = 15))
# 

save(list = c("slope_stanfit",
              "out_base",
              "stan_data",
              "jags_data",
              "route_map",
              "realized_strata_map",
              "firstYear",
              "sp_file",
              "species_f",
              "csv_files",
              "output_dir"),
     file = sp_file)

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



for(species in allspecies.eng){
  
  
  species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)
  
  sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")
  
  
  #sp_file <- paste0("output/",species,"Canadian_",firstYear,"_",lastYear,"_slope_route_iCAR2.RData")
  if(file.exists(sp_file)){

    load(sp_file)
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



beta_samples <- beta_samples %>% 
  ungroup() %>% 
  select(s,.draw,.value) %>% 
  rename(beta = .value)
  
alpha_samples <- alpha_samples %>% 
  ungroup() %>% 
  select(s,.draw,.value) %>% 
  rename(alpha = .value)

ab_samples = inner_join(beta_samples,alpha_samples)

ab_samples = left_join(ab_samples,sdnoise_samples,by = ".draw")

nyears = stan_data$nyears
fixedyear = stan_data$fixedyear
YEARS = c(min(jags_data$r_year):max(jags_data$r_year))

if(length(YEARS) != nyears){stop("years don't match YEARS =",length(YEARS),"nyears =",nyears)}

ind_fxn = function(a,b,sdn,y,fy){
  i = exp(a + b*(y-fy) + (0.5*(sdn^2)))
  return(i)
}

i_samples = NULL
for(yr in 1:nyears){
  i_t = ab_samples %>% 
    mutate(i = ind_fxn(alpha,beta,sdnoise,yr,fixedyear),
           y = yr,
           year = YEARS[yr])
  i_samples <- bind_rows(i_samples,i_t)
}

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



# add mapping of trends ---------------------------------------------------


breaks <- c(-0.07, -0.04, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.04, 0.07)
labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
labls = paste0(labls, " slope")
map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                 "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
names(map_palette) <- labls

route_map_out$Tplot <- cut(route_map_out$b,breaks = c(-Inf, breaks, Inf),labels = labls)

tmap = ggplot(route_map_out)+
  geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(aes(colour = Tplot,size = abund))+
  scale_size_continuous(range = c(0.5,3))+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0("slope\n",firstYear,"-",lastYear))+
  labs(title = paste(species,"trend-slopes by route (size = mean abundance)"))

maps[[jj]] <- tmap

tmap2 = ggplot(route_map_out)+
  geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(aes(colour = Tplot,size = abund))+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0("slope\n",firstYear,"-",lastYear))+
  theme(legend.position = "none")+
  labs(title = paste(species))

maps2[[jj]] <- tmap2


route_map_out <- route_map_out %>% 
  mutate(h_ci = (uci-lci)/2)
tmap3 = ggplot(route_map_out)+
  geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(aes(colour = Tplot,size = 1/h_ci))+
  scale_size_continuous(range = c(0.5,3))+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0("slope\n",firstYear,"-",lastYear))+
  labs(title = paste(species))

maps3[[jj]] <- tmap3



route_map_out_space$Tplot <- cut(route_map_out_space$b,breaks = c(-Inf, breaks, Inf),labels = labls)

tmap_space = ggplot(route_map_out_space)+
  geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(aes(colour = Tplot,size = abund))+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0("slope\n",firstYear,"-",lastYear))+
  theme(legend.position = "none")+
  labs(title = paste("spatial component"))
maps_space[[jj]] <- tmap_space




route_map_out_rand$Tplot <- cut(route_map_out_rand$b,breaks = c(-Inf, breaks, Inf),labels = labls)

tmap_rand = ggplot(route_map_out_rand)+
  geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(aes(colour = Tplot,size = abund))+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0("slope\n",firstYear,"-",lastYear))+
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

write.csv(trends_out,
          file = paste0("output/combined_",firstYear,"_",lastYear,"_",scope,"_trends_and_intercepts2.csv"))



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



write.csv(trends_out_space,
          file = paste0("output/combined_",firstYear,"_",lastYear,"_",scope,"_spatial_trends_and_intercepts2.csv"))
write.csv(trends_out_rand,
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

## building a Stan version of the bbsBayes models

library(bbsBayes)
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE,javascript = FALSE)
library(shinystan)
library(sf)
library(spdep)
library(doParallel)
library(foreach)
# library(ggforce)
# library(tidybayes)
source("functions/mungeCARdata4stan.R")
source("functions/prepare-jags-data-alt.R") ## small alteration of the bbsBayes function
## changes captured in a commit on Nov 20, 2020

pkgs = c("bbsBayes","tidyverse",
         "doParallel","foreach",
         "rstan","sf","spdep")

# load and stratify CASW data ---------------------------------------------
#species = "Pacific Wren"
#species = "Barred Owl"
strat = "bbs_usgs"
model = "slope"

strat_data = stratify(by = strat)

firstYear = 1995
lastYear = 2019

species_list = c("Bobolink",
                 "Eastern Meadowlark",
                 "Chimney Swift",
                 "Baird's Sparrow",
                 "Sprague's Pipit",
                 "Common Nighthawk",
                 "Canada Warbler",
                 "Blackpoll Warbler",
                 "Bank Swallow",
                 "Purple Martin",
                 "Barn Swallow",
                 "Brown-headed Cowbird",
                 "Vesper Sparrow",
                 "Red-winged Blackbird",
                 "Eastern Kingbird",
                 "Killdeer",
                 "Northern Harrier",
                 "Savannah Sparrow",
                 "Western Meadowlark",
                 "American Kestrel",
                 "Bobolink",
                 "Horned Lark",
                 "Grasshopper Sparrow",
                 "Baird's Sparrow",
                 "Eastern Meadowlark",
                 "Sprague's Pipit",
                 "Chestnut-collared Longspur",
                 "Golden-winged Warbler")


# optional parallel setup ----------------------------------------------------------
# n_cores <- 5
# cluster <- makeCluster(n_cores,type = "PSOCK")
# registerDoParallel(cluster)

# nspecies <- length(allspecies.eng)


# allsum <- foreach(species = species_list,
#                   .packages = pkgs,
#                   .inorder = FALSE,
#                   .errorhandling = "pass") %dopar% {

                    # for(ssi in which(allspecies.eng %in% speciestemp2)){
                    
#for(species in species_list){
  
species = species_list[2]


# setting up covariate data -----------------------------------------------
  
  
#getting list of the routes in the data

  rt_list = unique(strat_data$route_strat[,c("strat_name","rt.uni")])
  
  
  load("covariate_data/compiled_footprint_data.RData")
  # 
  # #select scale for intercepts
  # # buffer_sizes
  # # [1] "400m" "1km"  "2km"  "3km"  "4km"  "5km" 
  sc = buffer_sizes[2] #1km
  # 
  # 
  # #select canadian footprint variable
  # # fp_components
  # # [1] "cumulative"                   "built"                        "crop"                         "dam_and_associated_reservoir" "forestry_harvest"             "mines"                       
  # # [7] "nav_water"                    "night_lights"                 "oil_gas"                      "pasture"                      "population_density"           "rail"                        
  # # [13] "roads"
  
  
  
  # slope covariates --------------------------------------------------------
  
  #slope covariates - change in footprint uisng global footprint layers 2000 - 2013
  vl_val = "mean" #alternate is "p_above" = proportion of landscape above a threshold FP value
  cls = paste(vl_val,sc,sep = "_")
  
  covs_sel <- fp_global_by_route[,c("rt.uni","year",cls)]
  covs_sel[,"y1"] <- covs_sel[,cls]
  covs_sel[,"x1"] <- covs_sel[,"year"]
  
  chng_fxn = function(x,y){
    s = as.numeric(coefficients(lm(y ~ x))["x"])
    s = s*diff(range(x))
    return(s)
  }
  
  # tmp = covs_sel[which(covs_sel$rt.uni == "11-13"),]
  
  fp_change <- covs_sel %>% filter(!is.na(x1) & !is.na(y1)) %>% 
    group_by(rt.uni) %>% 
    summarise(Global_FP_Change = chng_fxn(x1,y1))
  
 
  # intercept covariates ----------------------------------------------------
  fp_covi = fp_components[c(1,10)] #cumulative and pasture
  covi_head = paste(fp_covi,sc,"mean",sep = "_")
  
  fpi = fp_can_by_route[,c("rt.uni",covi_head)]
  
  all_covs1 = inner_join(fpi,fp_change)
  
  all_covs = inner_join(rt_list,all_covs1,by = c("rt.uni"))
 
  
  ### identifying which covariates are for slope and intercepts
  cov_names_slopes = c(covi_head[1],"Global_FP_Change")
  
  cov_names_intercepts = covi_head
  
# Modifying the object created by bbsBayes stratify() function
# that is removing data from routes with no covariates ----------------------------

  strat_data$route_strat <- strat_data$route_strat[which(strat_data$route_strat$rt.uni %in% all_covs$rt.uni),]
  strat_data$bird_strat <- strat_data$bird_strat[which(strat_data$bird_strat$rt.uni %in% all_covs$rt.uni),]
  
  
  jags_data = prepare_jags_data(strat_data = strat_data,
                             species_to_run = species,
                             model = model,
                             #n_knots = 10,
                             min_year = firstYear,
                             max_year = lastYear,
                             min_n_routes = 1) # 
### now just the BBS data with covariates remain.

 
 #  #### alternative approach that specifically removes the US data
 #  #### not necessary here, because only the Canadian data have covariate info, so US data already removed
 #  
 #   # removing the non-Canadian data ------------------------------------------
 #  names_strata <- get_composite_regions(strata_type = strat) ## bbsBayes function that gets a list of the strata names and their composite regions (provinces, BCRs, etc.)
 #  
 #  us_strata_remove <- names_strata[which(names_strata$national == "US"),"region"] # character vector of the strata names to remove from data
 #  
 #  
 # jags_data = prepare_jags_data(strat_data = strat_data,
 #                                species_to_run = species,
 #                                model = model,
 #                                #n_knots = 10,
 #                                min_year = firstYear,
 #                                max_year = lastYear,
 #                                min_n_routes = 1,
 #                                strata_rem = us_strata_remove) # this final argument removes all data from the US
 #  ### now just hte Canadian data remain.
 #  
# spatial neighbourhood define --------------------------------------------
laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system

locat = system.file("maps",
                    package = "bbsBayes")
map.file = "BBS_USGS_strata"

strata_map = read_sf(dsn = locat,
                     layer = map.file)
strata_map = st_transform(strata_map,crs = laea)

real_strata_map = filter(strata_map,ST_12 %in% unique(jags_data$strat_name))

strata_list <- data.frame(ST_12 = unique(jags_data$strat_name),
                          strat = unique(jags_data$strat))


real_strata_map <- inner_join(real_strata_map,strata_list, by = "ST_12")
strata_bounds <- st_union(real_strata_map)
strata_bounds_buf = st_buffer(strata_bounds,dist = 100000)


jags_data[["routeF"]] <- as.integer(factor((jags_data$route)))

route_map = unique(data.frame(route = jags_data$route,
                              routeF = jags_data$routeF,
                              strat = jags_data$strat_name,
                              Latitude = jags_data$Latitude,
                              Longitude = jags_data$Longitude))

all_covs <- inner_join(route_map,all_covs,by = c("route" = "rt.uni")) %>% 
  arrange(routeF)


if(nrow(all_covs) != max(route_map$routeF)){stop("Error - missing bird or covariate data")}


cov_names = names(all_covs)[7:ncol(all_covs)]

### scales the predictors, 

for(j in 1:length(cov_names)){
  j1 = cov_names[j]
  all_covs[,paste0("sc_",j1)] = scale(all_covs[,j1]) ## scales the linear predictor
  all_covs[,paste0("sc_",j1,"_2")] = (all_covs[,paste0("sc_",j1)]^2)-mean(unlist(all_covs[,paste0("sc_",j1)]^2)) ## creates the polynomial, then centers
}

# merge the covariates with BBS data --------------------------------------

beta_covs <- as.matrix((all_covs[,paste0("sc_",cov_names_slopes)]))

beta_covs2 <- as.matrix(all_covs[,paste0("sc_",cov_names_slopes,"_2")])
n_c_beta <- length(cov_names_slopes)




alpha_covs <- as.matrix((all_covs[,paste0("sc_",cov_names_intercepts)]))

alpha_covs2 <- as.matrix(all_covs[,paste0("sc_",cov_names_intercepts,"_2")])
n_c_alpha <- length(cov_names_intercepts)





# reconcile duplicate spatial locations -----------------------------------
dups = which(duplicated(route_map[,c("Latitude","Longitude")]))
while(length(dups) > 0){
  route_map[dups,"Latitude"] <- route_map[dups,"Latitude"]+0.01 #=0.01 decimal degrees ~ 1km
  route_map[dups,"Longitude"] <- route_map[dups,"Longitude"]+0.01 #=0.01 decimal degrees ~ 1km
  dups = which(duplicated(route_map[,c("Latitude","Longitude")]))
  
}
#dups = which(duplicated(route_map[,c("Latitude","Longitude")]))


route_map = st_as_sf(route_map,coords = c("Longitude","Latitude"))
st_crs(route_map) <- 4269 #NAD83 commonly used by US federal agencies
#load strata map


route_map = st_transform(route_map,crs = laea)

# generate neighbourhoods -------------------------------------------------

# coords = st_coordinates(route_map)


# nb_dbpoly <- spdep::poly2nb(real_strata_map,row.names = real_strata_map$strat)
# 
# plot(real_strata_map,max.plot = 1,reset = F)
# plot(nb_dbpoly,coords,add = T,col = "red")


# Voronoi polygons from route locations -----------------------------------
box <- st_as_sfc(st_bbox(route_map))

v <- st_cast(st_voronoi(st_union(route_map), envelope = box))

vint = st_sf(st_cast(st_intersection(v,strata_bounds_buf),"POLYGON"))
vintj = st_join(vint,route_map,join = st_contains)
vintj = arrange(vintj,routeF)

nb_db = poly2nb(vintj,row.names = vintj$route,queen = FALSE)


### currently using 2 nearest neighbours to define the spatial relationships
## many regions may  have > 2 neighbours because of the symmetry condition
# nb_db <- spdep::knn2nb(spdep::knearneigh(coords,k = 4),row.names = route_map$route,sym = TRUE)
cc = st_coordinates(st_centroid(vintj))
#

ggp = ggplot(data = route_map)+
  geom_sf(data = vintj,alpha = 0.3)+
  geom_sf(aes(col = strat))+
  geom_sf_text(aes(label = routeF),size = 3,alpha = 0.3)+
  theme(legend.position = "none")
pdf(file = paste0("output/",species,"Canadian route maps ",firstYear," ",lastYear,".pdf"))
plot(nb_db,cc,col = "red")
print(ggp)
dev.off()


# wca = which(grepl(route_map$strat,pattern = "-CA-",fixed = T))
# wak = which(grepl(route_map$strat,pattern = "-AK-",fixed = T))
# 
# nb2[[wak]]


nb_info = spdep::nb2WB(nb_db)

### re-arrange GEOBUGS formated nb_info into appropriate format for Stan model
car_stan_dat <- mungeCARdata4stan(adjBUGS = nb_info$adj,
                                  numBUGS = nb_info$num)




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


stan_data[["beta_covs"]] <- beta_covs
#stan_data[["beta_covs2"]] <- beta_covs2
stan_data[["n_c_beta"]] <- n_c_beta


stan_data[["alpha_covs"]] <- alpha_covs
#stan_data[["alpha_covs2"]] <- alpha_covs2
stan_data[["n_c_alpha"]] <- n_c_alpha


if(car_stan_dat$N != stan_data[["nroutes"]]){stop("Some routes are missing from adjacency matrix")}

mod.file = "models/slope_iCAR_route_covariates.stan"




# model fitting -----------------------------------------------------------


parms = c("sdnoise",
          "sdobs",
          "sdbeta",
          "alpha",
          "sdalpha",
          "BETA",
          "ALPHA",
          "beta",
          "eta",
          "log_lik",
          "c_beta",
          #"c_beta2",
          "sum_beta")

## compile model
slope_model = stan_model(file=mod.file)

print(species)
## run sampler on model, data
slope_stanfit <- sampling(slope_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=3, iter=1300,
                               warmup=1000,
                               cores = 3,
                               pars = parms,
                               control = list(adapt_delta = 0.99,
                                              max_treedepth = 15))


save(list = c("slope_stanfit","stan_data","jags_data","vintj","route_map","real_strata_map"),
     file = paste0("output/",species,"Canadian_",firstYear,"_",lastYear,"_covariate_route_iCAR.RData"))

#}


load(paste0("output/",species,"Canadian_",firstYear,"_",lastYear,"_covariate_route_iCAR.RData"))


#stopCluster(cl = cluster)








# post loop analysis ------------------------------------------------------


# 
 launch_shinystan(slope_stanfit) 
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


# plotting and trend output -----------------------------------------------

library(tidybayes)


for(species in species_list){
  
  sp_file <- paste0("output/",species,"Canadian_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")
  if(file.exists(sp_file)){

    load(sp_file)
    ### may be removed after re-running
    
    laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system
    
    locat = system.file("maps",
                        package = "bbsBayes")
    map.file = "BBS_USGS_strata"
    
    strata_map = read_sf(dsn = locat,
                         layer = map.file)
    strata_map = st_transform(strata_map,crs = laea)
    
    real_strata_map = filter(strata_map,ST_12 %in% unique(jags_data$strat_name))
    
    strata_list <- data.frame(ST_12 = unique(jags_data$strat_name),
                              strat = unique(jags_data$strat))
    
    
    real_strata_map <- inner_join(real_strata_map,strata_list, by = "ST_12")
    
    
    ####
# add trend and abundance ----------------------------------------

beta_samples = gather_draws(slope_stanfit,beta[s])

slopes = beta_samples %>% group_by(s) %>% 
  summarise(b = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975),
            sd = sd(.value),
            prec = 1/var(.value),
            .groups = "keep")

alpha_samples = gather_draws(slope_stanfit,alpha[s])
interc = alpha_samples %>% group_by(s) %>% 
  summarise(abund = mean(exp(.value)),
            lci_i = quantile(exp(.value),0.025),
            uci_i = quantile(exp(.value),0.975),
            sd_i = sd(exp(.value)),
            prec_i = 1/var(.value),
            .groups = "keep")

#plot(log(interc$i),slopes$b)
slops_int = inner_join(slopes,interc,by = "s")
slops_int$routeF = slops_int$s




# connect trends to original route names ----------------------------------

route_map_out = left_join(route_map,slops_int,by = "routeF")

# add mapping of trends ---------------------------------------------------


breaks <- c(-0.07, -0.04, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.04, 0.07)
labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
labls = paste0(labls, " slope")
route_map_out$Tplot <- cut(route_map_out$b,breaks = c(-Inf, breaks, Inf),labels = labls)
map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                 "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
names(map_palette) <- labls


tmap = ggplot(route_map_out)+
  geom_sf(data = real_strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(aes(colour = Tplot,size = abund))+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0("slope\n",firstYear,"-",lastYear))+
  labs(title = paste(species,"trend-slopes by route (size = mean abundance)"))
# Send to Courtney --------------------------------------------------------


pdf(file = paste0("figures/",species,firstYear,"_",lastYear,"_Canadian_trend_map_route.pdf"),
    height = 8.5,
    width = 11)

print(tmap)

dev.off()

write.csv(route_map_out,
          file = paste0("output/",species," ",firstYear," ",lastYear,"_Canadian_trends_and_intercepts.csv"))


  }
}


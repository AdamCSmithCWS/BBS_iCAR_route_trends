## building a Stan version of the bbsBayes models

library(bbsBayes)
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE, javascript = FALSE)
library(shinystan)
library(sf)
library(spdep)
library(doParallel)
library(foreach)
 library(ggforce)
library(tidybayes)
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

firstYear = 2004
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

allspecies.eng = strat_data$species_strat$english

species_list = allspecies.eng[-which(allspecies.eng %in% species_list)]



# removing the non-Canadian data ------------------------------------------

names_strata <- get_composite_regions(strata_type = strat) ## bbsBayes function that gets a list of the strata names and their composite regions (provinces, BCRs, etc.)

us_strata_remove <- names_strata[which(names_strata$national == "US"),"region"] # character vector of the strata names to remove from data



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
                    
for(species in rev(allspecies.eng)){
  
  sp_file <- paste0("output/",species,"Canadian_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")
  if(file.exists(sp_file)){next}
    
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

if(min(nb_info$num) == 0){next}
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


if(car_stan_dat$N != stan_data[["nroutes"]]){stop("Some routes are missing from adjacency matrix")}

mod.file = "models/slope_iCAR_route.stan"

parms = c("sdnoise",
          "sdobs",
          "sdbeta",
          "alpha",
          "sdalpha",
          "BETA",
          "ALPHA",
          "beta",
          "eta",
          "log_lik")

## compile model
slope_model = stan_model(file=mod.file)

print(paste(firstYear,species))
## run sampler on model, data
slope_stanfit <- sampling(slope_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=4, iter=900,
                               warmup=600,
                               cores = 4,
                               pars = parms,
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 15))


save(list = c("slope_stanfit","stan_data","jags_data","vintj","route_map","real_strata_map"),
     file = paste0("output/",species,"Canadian_",firstYear,"_",lastYear,"_slope_route_iCAR.RData"))

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

maps = vector(mode = "list",length = 300)
jj <- 0
trends_out <- NULL

for(species in allspecies.eng){
  
  sp_file <- paste0("output/",species,"Canadian_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")
  if(file.exists(sp_file)){

    load(sp_file)
    ### may be removed after re-running
    jj <- jj+1
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




# Route-level trajectories ------------------------------------------------
if(route_trajectories){
sdnoise_samples = gather_draws(slope_stanfit,sdnoise)%>% 
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
            lci_i = quantile(i,0.025),
            uci_i = quantile(i,0.975),
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
pdf(paste0("trajectories/",species,"_route_trajectories.pdf"),
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




maps[[jj]] <- tmap


# write.csv(route_map_out,
#           file = paste0("output/",species," ",firstYear," ",lastYear,"_Canadian_trends_and_intercepts.csv"))


  }
}



pdf(file = paste0("figures/Combined_",firstYear,"_",lastYear,"_Canadian_trend_map_route.pdf"),
    height = 8.5,
    width = 11)
for(j in 1:length(maps)){
  if(!is.null(maps[[j]])){print(maps[[j]])}
}
dev.off()

write.csv(trends_out,
          file = paste0("output/combined_",firstYear,"_",lastYear,"_Canadian_trends_and_intercepts.csv"))




## Fitting the BYM model to 1995 - 2021 BBS data
## script currently written to fit the model then save the Stan output to a directory
## 
#setwd("C:/GitHub/BBS_iCAR_route_trends")
library(bbsBayes)
library(tidyverse)
library(cmdstanr)
source("functions/neighbours_define.R") ## function to define neighbourhood relationships
source("functions/prepare-data-alt.R") ## small alteration of the bbsBayes function
## above source() call over-writes the bbsBayes prepare_data() function
source("functions/get_basemap_function.R") ## loads one of the bbsBayes strata maps
source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output


# load and stratify CASW data ---------------------------------------------
#species = "Pacific Wren"
strat = "bbs_usgs"
model = "slope"


firstYear = 1995
lastYear = 2021
output_dir <- "output"

## this list should include all of the species that we're interested in for the grasslands project
species_list = "Dickcissel"


# these two can be set to FALSE if it makes sense just to fit the models and store the model output from Stan
# e.g., if it's easier to run the models and then to produce maps and trends after the fact
produce_trends <- TRUE
produce_maps <- TRUE # this is only relevant if produce_trends == TRUE

# I've got this running as a simple species loop
# it would be more efficient to run it in parallel using the foreach and parallel packages, but I can't seem to get Stan to work using these parallel options
for(species in species_list){
  
species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)

spp <- "_BYM_"

out_base <- paste0(species_f,spp,firstYear,"_",lastYear)


# SPECIES LOOP ------------------------------------------------------------




sp_data_file <- paste0(output_dir,"/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")


strat_data <- bbsBayes::stratify(by = strat)

### this is the alternate prepare data function
jags_data = prepare_data(strat_data = strat_data,
                              species_to_run = species,
                              model = model,
                              #n_knots = 10,
                              min_year = firstYear,
                              max_year = lastYear,
                              min_n_routes = 1)# spatial neighbourhood define --------------------------------------------

# strata map of one of the bbsBayes base maps
# helps group and set boundaries for the route-level neighbours,
## NOT directly used in the model
strata_map  <- get_basemap(strata_type = strat,
                           transform_laea = TRUE,
                           append_area_weights = FALSE)

strata_list <- data.frame(ST_12 = unique(jags_data$strat_name),
                          strat = unique(jags_data$strat))

#drops the strata that are not included for this species
# and adds the strat indicator variable to link to model output
realized_strata_map <- filter(strata_map,ST_12 %in% unique(jags_data$strat_name)) %>% 
  inner_join(.,strata_list, by = "ST_12")

# Spatial boundaries set up --------------------

jags_data[["routeF"]] <- as.integer(factor((jags_data$route)))

#create a data frame of each unique route in the species-specific dataset
route_map = unique(data.frame(route = jags_data$route,
                              routeF = jags_data$routeF,
                              strat = jags_data$strat_name,
                              Latitude = jags_data$Latitude,
                              Longitude = jags_data$Longitude))


# reconcile duplicate spatial locations -----------------------------------
# adhoc way of separating different routes with the same starting coordinates
# this shifts the starting coordinates of teh duplicates by ~1.5km to the North East 
# ensures that the duplicates have a unique spatial location, but remain very close to
# their original location and retain reasonable neighbourhood relationships
# these duplicates happen when a "new" route is established because some large proportion
# of the end of a route is changed, but the start-point remains the same
dups = which(duplicated(route_map[,c("Latitude","Longitude")]))
while(length(dups) > 0){
  route_map[dups,"Latitude"] <- route_map[dups,"Latitude"]+0.01 #=0.01 decimal degrees ~ 1km
  route_map[dups,"Longitude"] <- route_map[dups,"Longitude"]+0.01 #=0.01 decimal degrees ~ 1km
  dups = which(duplicated(route_map[,c("Latitude","Longitude")]))
  
}
dups = which(duplicated(route_map[,c("Latitude","Longitude")])) 
if(length(dups) > 0){stop(paste(spec,"At least one duplicate route remains"))}


#create spatial object from route_map dataframe
route_map = st_as_sf(route_map,coords = c("Longitude","Latitude"))
st_crs(route_map) <- 4269 #NAD83 commonly used by US federal agencies

#reconcile the projections of routes and base bbs strata
route_map = st_transform(route_map,crs = st_crs(realized_strata_map))


## custom function that returns the adjacency data necessary for the stan model
## also exports maps and saved data objects to plot_dir
if(!dir.exists("data/maps/")){
  if(!dir.exists("data/")){
    dir.create("data/")
  }
  dir.create("data/maps/")
}
car_stan_dat <- neighbours_define(real_strata_map = route_map,
                                  #strat_link_fill = 100000,
                                  plot_neighbours = TRUE,
                                  species = species,
                                  plot_dir = "data/maps/",
                                  plot_file = paste0("_",firstYear,"_",lastYear,"_route_maps"),
                                  save_plot_data = TRUE,
                                  voronoi = TRUE,
                                  strat_indicator = "routeF",
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
     file = sp_data_file)






  mod.file = "models/slope_BYM_route_NB.stan"
  
  init_def <- function(){ list(alpha_raw = rnorm(stan_data$nroutes,0,0.1),
                               sdalpha = runif(1,0.01,0.1),
                               ALPHA = 0,
                               BETA = 0,
                               eta = 0,
                               obs_raw = rnorm(stan_data$nobservers,0,0.1),
                               sdnoise = 0.2,
                               sdobs = 0.1,
                               sdbeta_rand = runif(1,0.01,0.1),
                               beta_raw_rand = rnorm(stan_data$nroutes,0,0.01),
                               sdbeta_space = runif(1,0.01,0.1),
                               beta_raw_space = rnorm(stan_data$nroutes,0,0.01))} 



slope_model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}
stanfit <- slope_model$sample(
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

summ <- stanfit$summary()

save(list = c("stanfit","summ"),
     file = paste0(output_dir,"/",out_base,"_stanfit.RData"))



# extract trends and abundances -------------------------------------------
if(produce_trends){
routes_df <- data.frame(routeF = jags_data$routeF,
                        route = jags_data$route,
                        year = jags_data$r_year,
                        obs = jags_data$obser,
                        stratum = jags_data$strat_name) %>% 
  group_by(route,routeF,obs,stratum) %>%
  summarise(n_yr_obs = n(),
            .groups = "drop") %>% 
  group_by(route,routeF,stratum) %>% 
  summarise(n_obs = n(),
            mean_y_obs = mean(n_yr_obs),
            max_y_obs = max(n_yr_obs),
            .groups = "drop") %>% 
  distinct()

# function to turn a log-scale slope parameter into a trend estimate in %/year
tr_f <- function(x){
  t <- (exp(x)-1)*100
}


# route-level trends and abundances ------------------------------------------


abunds <- posterior_samples(fit = stanfit,
                                  parm = "alpha",
                                  dims = "routeF") %>% 
  posterior_sums(.,
                 dims = "routeF")%>% 
  left_join(.,routes_df,by = "routeF") %>% 
  mutate(abund = exp(mean),
         abund_lci = exp(lci),
         abund_uci = tr_f(uci),
         abund_se = tr_f(sd),
         abund_Wci = abund_uci-abund_lci)%>% 
  select(route,abund,abund_lci,abund_uci,abund_se,abund_Wci)

trends <- posterior_samples(fit = stanfit,
                             parm = "beta",
                             dims = "routeF") %>%
  posterior_sums(.,
                 dims = "routeF")%>% 
  left_join(.,routes_df,by = "routeF") %>% 
  mutate(trend = tr_f(mean),
         trend_lci = tr_f(lci),
         trend_uci = tr_f(uci),
         trend_se = tr_f(sd),
         trend_Wci = trend_uci-trend_lci,
         log_scale_slope = mean,
         log_scale_slope_lci = lci,
         log_scale_slope_uci = uci,
         log_scale_slope_se = sd,
         log_scale_slope_Wci = log_scale_slope_uci-log_scale_slope_lci) %>% 
  select(route,trend,trend_lci,trend_uci,trend_se,trend_Wci,
         log_scale_slope,log_scale_slope_lci,log_scale_slope_uci,log_scale_slope_se,log_scale_slope_Wci)

est_table <- inner_join(trends,
                        abunds,
                        by = "route") %>% 
  mutate(species = species,
         trend_time = paste(firstYear,lastYear,sep = "-"))

if(!dir.exists("trends")){
  dir.create("trends")
}
#writes a csv file for each species
write.csv(est_table,
          paste0("trends/",species_f,"_trends_",firstYear,"-",lastYear,".csv"),
          row.names = FALSE)
#also appends species results to an all trends file
write.csv(est_table,
          paste0("trends/","All_trends_",firstYear,"-",lastYear,".csv"),
          row.names = FALSE,
          append = TRUE)
# if it makes sense to output the pdf maps of species trends
if(produce_maps){
  
route_bounds <- st_union(route_map) #union to provide a simple border of the realised strata
bb = st_bbox(route_bounds)
xlms = as.numeric(c(bb$xmin,bb$xmax))
ylms = as.numeric(c(bb$ymin,bb$ymax))


# merging route map with the estimates for this species -------------------
trend_plot_map <- route_map %>% 
  left_join(.,est_table,by = "route") 


breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
lgnd_head <- "Mean Trend\n"
trend_title <- "Mean Trend"
labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
labls = paste0(labls, " %/year")
trend_plot_map$Tplot <- cut(trend_plot_map$trend,breaks = c(-Inf, breaks, Inf),labels = labls)


map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                 "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
names(map_palette) <- labls


tmap = ggplot(trend_plot_map)+
  #geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(aes(colour = Tplot,size = abund))+
  scale_size_continuous(range = c(0.05,2),
                        name = "Mean Count")+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0(lgnd_head,firstYear,"-",lastYear))+
  coord_sf(xlim = xlms,ylim = ylms)

if(!file.exists("Figures")){
  file.create("Figures")
}
pdf(file = paste0("Figures/",species_f,"_",firstYear,"-",lastYear,"_trend_map.pdf"),
    width = 10,
    height = 8)
print(tmap)
dev.off()

}#end if produce_maps 

}# end if produce_trends

}# end of species loop



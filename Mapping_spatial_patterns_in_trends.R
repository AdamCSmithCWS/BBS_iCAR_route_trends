
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
names(trends_out2) <- c(names(trends_out2)[-1],"geometry2")
names(trends_out_space2) <- c(names(trends_out_space2)[-1],"geometry2")



sps_all <- unique(trends_out_space2$english_name)

socb <- read.csv("data/SOCB data supplement.csv")

groups = c("Grassland.birds",
           "Forest.birds",
           "Other.birds",
           "Aerial.insectivores",
           "suburban",
           "other.wetland.birds")



# group loop --------------------------------------------------------------

g_sel = groups[1]

for(g_sel in groups){

# species selection of trend data -----------------------------------------


socb_c = which(grepl(names(socb),pattern = g_sel))

sp_sel1 = socb[which(socb[,socb_c] == "Included in group"),"species"]
sp_sel1[which(sp_sel1 == "Le Conte's Sparrow")] <- "LeConte's Sparrow"
sp_sel2 = sp_sel1[which(sp_sel1 %in% sps_all)]



dat = trends_out_space2 %>% filter(english_name %in% sp_sel2)

#Trend = as.numeric(scale(Trend,scale = FALSE)),  


dat_plot <- dat %>% group_by(english_name) %>% 
  mutate(Abund = as.numeric(scale(Mean_abundance,scale = TRUE))) %>% 
  group_by(BBS_route,BBS_stratum) %>% 
  summarise(mean_Trend = mean(Trend),
            mean_Abund = mean(Abund),
            n_species = n())



# generate route map ------------------------------------------------------
strata_map  <- bbsBayes::load_map(stratify_by = strat)


route_map <- unique(strat_data$route_strat[,c("rt.uni","Latitude","Longitude","strat_name")])

route_map = st_as_sf(route_map,coords = c("Longitude","Latitude"))
st_crs(route_map) <- 4269 #NAD83 commonly used by US federal agencies
#load strata map

route_map = st_transform(route_map,crs = st_crs(strata_map))


route_map_out = inner_join(route_map,dat_plot,by = c("rt.uni" = "BBS_route",
                                            "strat_name" = "BBS_stratum"))



strata_bounds <- st_union(route_map_out) #union to provide a simple border of the realised strata
bb = st_bbox(strata_bounds)
xlms = as.numeric(c(bb$xmin,bb$xmax))
ylms = as.numeric(c(bb$ymin,bb$ymax))


# MAPPING -----------------------------------------------------------------

plot_trend <- TRUE #set to false to plot the abundance
if(plot_trend){
  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  lgnd_head <- "Mean Trend\n"
  trend_title <- "Mean Trend"
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %/year")
  route_map_out$Tplot <- cut(route_map_out$mean_Trend,breaks = c(-Inf, breaks, Inf),labels = labls)

  
  
}else{
  breaks <- c(-2, -1, -0.5, -0.1, -0.05, 0.05, 0.1, 0.5, 1, 2)
  lgnd_head <- "Mean Scaled Abundance\n"
  trend_title <- "Mean Scaled Abundance"
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " Abund")
  route_map_out$Tplot <- cut(route_map_out$mean_Abund,breaks = c(-Inf, breaks, Inf),labels = labls)

}
map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                 "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
names(map_palette) <- labls


tmap = ggplot(route_map_out)+
  #geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(aes(colour = Tplot,size = n_species))+
  scale_size_continuous(range = c(0.1,3),
                        name = "Number of species")+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0(lgnd_head,firstYear,"-",lastYear))+
  coord_sf(xlim = xlms,ylim = ylms)+
  labs(title = paste("DRAFT ",g_sel,trend_title,"by BBS route"),
       subtitle = "Route-level trends from a spatial iCAR model, using Stan")



png(filename = paste0("Figures/images/",g_sel,"_Mean_Trends_",firstYear,".png"),
    res = 600,
    width = 20,
    height = 15,
    units = "cm")
print(tmap)
dev.off()



}





# just a standard route-map -----------------------------------------------
strata_map  <- bbsBayes::load_map(stratify_by = strat)


route_map <- unique(strat_data$route_strat[,c("rt.uni","Latitude","Longitude","strat_name")])

route_map = st_as_sf(route_map,coords = c("Longitude","Latitude"))
st_crs(route_map) <- 4269 #NAD83 commonly used by US federal agencies
#load strata map

route_map = st_transform(route_map,crs = st_crs(strata_map))

strata_bounds <- st_union(route_map) #union to provide a simple border of the realised strata
bb = st_bbox(strata_bounds)
xlms = as.numeric(c(bb$xmin,bb$xmax))
ylms = as.numeric(c(bb$ymin,bb$ymax))


# MAPPING -----------------------------------------------------------------



tmap = ggplot(route_map)+
  #geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(size = 0.1)+
  coord_sf(xlim = xlms,ylim = ylms)



png(filename = "all_bbs_routes.png",
    res = 600,
    width = 20,
    height = 15,
    units = "cm")
print(tmap)
dev.off()








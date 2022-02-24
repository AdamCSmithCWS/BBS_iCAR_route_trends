
# spatial patterns in mean trends among species ---------------------------


library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)
library(sf)
library(spdep)
library(ggforce)
#library(tidybayes)
#source("functions/mungeCARdata4stan.R")
source("functions/neighbours_define.R") ## function to define neighbourhood relationships
source("functions/prepare-jags-data-alt.R") ## small alteration of the bbsBayes function
source("functions/get_basemap_function.R") ## loads one of the bbsBayes strata maps
source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output
## changes captured in a commit on Nov 20, 2020


strat = "bbs_usgs"
model = "slope"


firstYear = 2004
lastYear = 2019

scope = "RangeWide"
check_conv = FALSE #set to TRUE to run convergence summary

species = "Blue-headed Vireo"
species = "Dickcissel"
species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)




# SPECIES LOOP ------------------------------------------------------------

output_dir <- "output"
sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")

load(sp_file) 

output_dir <- "output"

out_base_ns <- paste0(species_f,"_Non_spatial_",firstYear,"_",lastYear)

csv_files_ns <- paste0(output_dir,"/",out_base_ns,"-",1:3,".csv")

stanfit_ns <- as_cmdstan_fit(files = csv_files_ns)

if(check_conv){
  
stanf_df_ns <- stanfit_ns$draws(format = "df")


conv_summ <- summarise_draws(stanf_df_ns) %>% 
  mutate(species = species,
         model = out_base_ns)

}
out_base <- paste0(species_f,"_RangeWide_",firstYear,"_",lastYear)

csv_files <- paste0(output_dir,"/",out_base,"-",1:3,".csv")

stanfit <- as_cmdstan_fit(files = csv_files)

if(check_conv){
  
stanf_df <- stanfit$draws(format = "df")


conv_summ_out <- summarise_draws(stanf_df_ns) %>% 
  mutate(species = species,
         model = out_base) %>% 
  bind_rows(conv_summ)


failed_rhat <- conv_summ_out %>% 
  filter(rhat > 1.05)

failed_ess_bulk <- conv_summ_out %>% 
  filter(ess_bulk < 100)


write.csv(conv_summ_out,file = paste0("trends_etc/conv_summ_",out_base,".csv"))
}
# extract trends and abundances -------------------------------------------

routes_df <- data.frame(routeF = jags_data$routeF,
                      route = jags_data$route,
                      year = jags_data$r_year,
                      obs = jags_data$obser) %>% 
  group_by(route,routeF,obs) %>%
  summarise(n_yr_obs = n(),
            .groups = "drop") %>% 
  group_by(route,routeF) %>% 
  summarise(n_obs = n(),
            mean_y_obs = mean(n_yr_obs),
            max_y_obs = max(n_yr_obs),
            .groups = "drop") %>% 
  distinct()

tr_f <- function(x){
  t <- (exp(x)-1)*100
}

trends <- posterior_samples(fit = stanfit,
                              parm = "beta",
                              dims = "routeF") %>%
  posterior_sums(.,
                 dims = "routeF")%>% 
  left_join(.,routes_df,by = "routeF") %>% 
  mutate(trend = tr_f(mean),
         trend_lci = tr_f(lci),
         trend_uci = tr_f(uci),
         trend_se = tr_f(sd))





trends_ns <- posterior_samples(fit = stanfit_ns,
                            parm = "beta",
                            dims = "routeF") %>%
  posterior_sums(.,
                 dims = "routeF")%>% 
  left_join(.,routes_df,by = "routeF") %>% 
  mutate(trend = tr_f(mean),
         trend_lci = tr_f(lci),
         trend_uci = tr_f(uci),
         trend_se = tr_f(sd))

# abundances --------------------------------------------------------------


abund <- posterior_samples(fit = stanfit,
                            parm = "alpha",
                            dims = "routeF") %>%
  posterior_sums(.,
                 dims = "routeF")%>% 
  left_join(.,routes_df,by = "routeF") %>% 
  mutate(abund = exp(mean),
         abund_lci = exp(lci),
         abund_uci = exp(uci),
         abund_se = exp(sd)) %>% 
  select(route,abund,abund_lci,abund_uci,abund_se)


abund_ns <- posterior_samples(fit = stanfit_ns,
                              parm = "alpha",
                              dims = "routeF") %>%
  posterior_sums(.,
                 dims = "routeF")%>% 
  left_join(.,routes_df,by = "routeF") %>% 
  mutate(abund = exp(mean),
         abund_lci = exp(lci),
         abund_uci = exp(uci),
         abund_se = exp(sd)) %>% 
  select(route,abund,abund_lci,abund_uci,abund_se)


abund_comp <- abund_ns %>% 
  rename(abund_ns = abund,
         abund_lci_ns = abund_lci,
         abund_uci_ns = abund_uci,
         abund_se_ns = abund_se) %>% 
  right_join(.,abund,by = "route") %>% 
  mutate(abund_dif = exp(log(abund)-log(abund_ns)),
         abund_se_dif = exp(log(abund_se) - log(abund_se_ns)))

# summarise the random component of the spatial beta --------
# use the random compnent as a scaling on the beta comparison
# large values of random == strong data that does not fit the spatial pattern

trends_rand <- posterior_samples(fit = stanfit,
                            parm = "beta_rand",
                            dims = "routeF") %>%
  posterior_sums(.,
                 dims = "routeF")%>% 
  left_join(.,routes_df,by = "routeF") %>% 
  mutate(trend = tr_f(mean),
         trend_lci = tr_f(lci),
         trend_uci = tr_f(uci),
         trend_se = tr_f(sd),
         abs_trend = abs(trend))%>% 
  select(route,trend,trend_lci,trend_uci,trend_se,abs_trend) %>% 
  rename_with(.,~paste0(.x,"_rand"),.cols = contains("trend")) 
 
trends_space <- posterior_samples(fit = stanfit,
                                 parm = "beta_space",
                                 dims = "routeF") %>%
  posterior_sums(.,
                 dims = "routeF")%>% 
  left_join(.,routes_df,by = "routeF") %>% 
  mutate(trend = tr_f(mean),
         trend_lci = tr_f(lci),
         trend_uci = tr_f(uci),
         trend_se = tr_f(sd),
         abs_trend = abs(trend))%>% 
  select(route,trend,trend_lci,trend_uci,trend_se,abs_trend) %>% 
  rename_with(.,~paste0(.x,"_space"),.cols = contains("trend")) 



trend_comp <- trends_ns %>% 
  select(route,trend,trend_lci,trend_uci,trend_se) %>% 
  rename_with(.,~paste0(.x,"_ns"),.cols = contains("trend")) %>% 
  right_join(.,trends,by = "route") %>% 
  right_join(.,trends_rand,by = "route") %>% 
  right_join(.,trends_space,by = "route") %>% 
  mutate(trend_dif = trend-trend_ns,
         trend_se_dif = trend_se - trend_se_ns)


tcplot <- ggplot(data = trend_comp,
                 aes(y = trend_ns,x = trend,colour = trend_se_ns))+
  geom_abline(slope = 1,intercept = 0)+
  geom_point(alpha = 0.5)+
  scale_colour_viridis_c(aesthetics = "colour",direction = 1)+
  theme_bw()

print(tcplot)
# generate route map ------------------------------------------------------
strata_map  <- bbsBayes::load_map(stratify_by = strat)
# 
# 
# route_map <- unique(strat_data$route_strat[,c("rt.uni","Latitude","Longitude","strat_name")])
# 
# route_map = st_as_sf(route_map,coords = c("Longitude","Latitude"))
# st_crs(route_map) <- 4269 #NAD83 commonly used by US federal agencies
# #load strata map
# 
# route_map = st_transform(route_map,crs = st_crs(strata_map))
# 
# 
# route_map_out = inner_join(route_map,dat_plot,by = c("rt.uni" = "BBS_route",
#                                             "strat_name" = "BBS_stratum"))
# 
# 

strata_bounds <- st_union(route_map) #union to provide a simple border of the realised strata
bb = st_bbox(strata_bounds)
xlms = as.numeric(c(bb$xmin,bb$xmax))
ylms = as.numeric(c(bb$ymin,bb$ymax))


trend_plot_map <- route_map %>% 
  left_join(.,trends,by = "route") %>% 
  left_join(.,abund,by = "route")

trend_plot_map_ns <- route_map %>% 
  left_join(.,trends_ns,by = "route") %>% 
  left_join(.,abund,by = "route")

# MAPPING -----------------------------------------------------------------

# plot_trend <- TRUE #set to false to plot the abundance
# if(plot_trend){
  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  lgnd_head <- "Mean Trend\n"
  trend_title <- "Mean Trend"
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %/year")
  trend_plot_map$Tplot <- cut(trend_plot_map$trend,breaks = c(-Inf, breaks, Inf),labels = labls)
  trend_plot_map_ns$Tplot <- cut(trend_plot_map_ns$trend,breaks = c(-Inf, breaks, Inf),labels = labls)
  
  
#   
# }else{
#   breaks <- c(-2, -1, -0.5, -0.1, -0.05, 0.05, 0.1, 0.5, 1, 2)
#   lgnd_head <- "Mean Scaled Abundance\n"
#   trend_title <- "Mean Scaled Abundance"
#   labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
#   labls = paste0(labls, " Abund")
#   route_map_out$Tplot <- cut(route_map_out$mean_Abund,breaks = c(-Inf, breaks, Inf),labels = labls)
# 
# }
map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                 "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
names(map_palette) <- labls


tmap = ggplot(trend_plot_map)+
  #geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(aes(colour = Tplot,size = abund))+
  scale_size_continuous(range = c(0.1,3),
                        name = "Mean Count")+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0(lgnd_head,firstYear,"-",lastYear))+
  coord_sf(xlim = xlms,ylim = ylms)+
  title("Spatial")

tmap_ns = ggplot(trend_plot_map_ns)+
  #geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(aes(colour = Tplot,size = abund))+
  scale_size_continuous(range = c(0.1,3),
                        name = "Mean Count")+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0(lgnd_head,firstYear,"-",lastYear))+
  coord_sf(xlim = xlms,ylim = ylms)+
  title("NonSpatial")

print(tmap + tmap_ns)+
  plot_layout(guides = "collect")

out = (tmap + tmap_ns)+plot_layout(nrow = 2,guides = "collect")

# png(filename = paste0("Figures/images/",g_sel,"_Mean_Trends_",firstYear,".png"),
#     res = 600,
#     width = 20,
#     height = 15,
#     units = "cm")
pdf(file = paste0("Figures/",species_f,out_base,"_vs_nonspatial.pdf"),
    height = 11,
    width = 8.5)
print(out)
dev.off()








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








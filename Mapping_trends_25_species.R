
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
check_conv = TRUE #set to TRUE to run convergence summary
conv_summary <- NULL



load("Data/sp_sel.RData")


for(species in sp_sel[1:5]){
  ##########
  ### Currently not working - needs to have full time-series models run
  ### with the same data-saving process as for example species
  ### remove the use of CV-based results, they don't include 2019
  #########
  
  species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)
  species_f <- gsub(species_f,pattern = "'",replacement = "",fixed = T)
  
  
output_dir <- "output"
sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")

# load(sp_file) 
load(paste0("data/",species_f,"CV_base_data.RData"))

output_dir <- "output"
# for(spp in c("_","_BYM_","_Non_spatial_")){
  out_base <- paste0(species_f,"_",firstYear,"_",lastYear)
#   

out_base_Non_spatial <- paste0(species_f,"_Non_spatial_",firstYear,"_",lastYear,"_CV")

csv_files_Non_spatial <- paste0(output_dir,"/",out_base_Non_spatial,"-",1:3,".csv")

stanfit_Non_spatial <- as_cmdstan_fit(files = csv_files_Non_spatial)

if(check_conv){
  
stanf_df_Non_spatial <- stanfit_Non_spatial$draws(format = "df")


conv_summ <- summarise_draws(stanf_df_Non_spatial) %>% 
  mutate(species = species,
         model = out_base_Non_spatial)
conv_summary <- bind_rows(conv_summary,conv_summ)

}
out_base_iCAR <- paste0(species_f,"_iCAR_",firstYear,"_",lastYear,"_CV")

csv_files_iCAR <- paste0(output_dir,"/",out_base_iCAR,"-",1:3,".csv")

stanfit_iCAR <- as_cmdstan_fit(files = csv_files_iCAR)

if(check_conv){
  
  stanf_df_iCAR <- stanfit_iCAR$draws(format = "df")
  
  
  conv_summ <- summarise_draws(stanf_df_iCAR) %>% 
    mutate(species = species,
           model = out_base_iCAR)
  conv_summary <- bind_rows(conv_summary,conv_summ)
  
}

out_base_BYM <- paste0(species_f,"_BYM_",firstYear,"_",lastYear,"_CV")

csv_files_BYM <- paste0(output_dir,"/",out_base_BYM,"-",1:3,".csv")

stanfit_BYM <- as_cmdstan_fit(files = csv_files_BYM)

if(check_conv){
  
  stanf_df_BYM <- stanfit_BYM$draws(format = "df")
  
  
  conv_summ <- summarise_draws(stanf_df__BYM) %>% 
    mutate(species = species,
           model = out_base_BYM)
  conv_summary <- bind_rows(conv_summary,conv_summ)
  



failed_rhat <- conv_summary %>% 
  filter(rhat > 1.05)

failed_ess_bulk <- conv_summary %>% 
  filter(ess_bulk < 100)


write.csv(conv_summary,file = paste0("trends_etc/conv_summ_",out_base,".csv"))
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

trend_long <- NULL
abund_long <- NULL
for(sppn in c("_iCAR","_BYM","_Non_spatial")){

  if(sppn == "_iCAR"){stanfit <- stanfit_iCAR}
  if(sppn == "_BYM"){stanfit <- stanfit__BYM}
  if(sppn == "_Non_spatial"){stanfit <- stanfit_Non_spatial}
  

# Trend gather ------------------------------------------------------------

  
trendst <- posterior_samples(fit = stanfit,
                              parm = "beta",
                              dims = "routeF") %>%
  posterior_sums(.,
                 dims = "routeF")%>% 
  left_join(.,routes_df,by = "routeF") %>% 
  mutate(trend = tr_f(mean),
         trend_lci = tr_f(lci),
         trend_uci = tr_f(uci),
         trend_se = tr_f(sd))%>% 
  select(route,trend,trend_lci,trend_uci,trend_se)

trend_ct <- trendst %>% 
  rename_with(.,~paste0(.x,sppn),.cols = contains("trend"))

trendst <- trendst %>% 
  mutate(version = sppn)

if(sppn == "_iCAR"){trend_comp <- trend_ct}else{
  trend_comp <- left_join(trend_comp,trend_ct,by = "route")
}
trend_long <- bind_rows(trend_long,trendst)



# abundance gather --------------------------------------------------------

abundst <- posterior_samples(fit = stanfit,
                             parm = "alpha",
                             dims = "routeF") %>%
  posterior_sums(.,
                 dims = "routeF")%>% 
  left_join(.,routes_df,by = "routeF") %>% 
  mutate(abund = exp(mean),
         abund_lci = exp(lci),
         abund_uci = exp(uci),
         abund_se = exp(sd))%>% 
  select(route,abund,abund_lci,abund_uci,abund_se)
abund_ct <- abundst %>% 
  rename_with(.,~paste0(.x,sppn),.cols = contains("abund"))

abundst <- abundst %>% 
  mutate(version = sppn)

if(sppn == "_iCAR"){abund_comp <- abund_ct}else{
  abund_comp <- left_join(abund_comp,abund_ct,by = "route")
}
abund_long <- bind_rows(abund_long,abundst)




}#end sppn


trends_rand <- posterior_samples(fit = stanfit__BYM,
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
 
trends_space <- posterior_samples(fit = stanfit__BYM,
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

save(list = c("trend_comp",
              "abund_comp",
              "trends_rand",
              "trends_space"),
     file = "Figures/example_trend_comparison_data.RData")


### compare the differences to the number of routes, and years with single obs etc.

tcplot <- ggplot(data = trend_comp,
                 aes(x = trend_iCAR,y = trend_Non_spatial,colour = trend_se_Non_spatial))+
  geom_abline(slope = 1,intercept = 0)+
  geom_point(alpha = 0.5)+
  scale_colour_viridis_c(aesthetics = "colour",direction = 1)+
  theme_bw()

tcplot2 <- ggplot(data = trend_comp,
                 aes(x = trend_iCAR,y = trend_BYM,colour = trend_se_BYM))+
  geom_abline(slope = 1,intercept = 0)+
  geom_point(alpha = 0.5)+
  scale_colour_viridis_c(aesthetics = "colour",direction = 1)+
  theme_bw()


print(tcplot + tcplot2)
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
  left_join(.,trend_long,by = "route") %>% 
  left_join(.,abund_long,by = c("route","version"))

save(list = c("trend_plot_map"),
     file = "Figures/example_trend_map_data.RData")
# trend_plot_map_ns <- route_map %>% 
#   left_join(.,trends_ns,by = "route") %>% 
#   left_join(.,abund,by = "route")

# MAPPING -----------------------------------------------------------------

# plot_trend <- TRUE #set to false to plot the abundance
# if(plot_trend){
  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  lgnd_head <- "Mean Trend\n"
  trend_title <- "Mean Trend"
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %/year")
  trend_plot_map$Tplot <- cut(trend_plot_map$trend,breaks = c(-Inf, breaks, Inf),labels = labls)

  
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
  scale_size_continuous(range = c(0.05,2),
                        name = "Mean Count")+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0(lgnd_head,firstYear,"-",lastYear))+
  coord_sf(xlim = xlms,ylim = ylms)+
  facet_wrap(vars(version),nrow = 3)

# png(filename = paste0("Figures/images/",g_sel,"_Mean_Trends_",firstYear,".png"),
#     res = 600,
#     width = 20,
#     height = 15,
#     units = "cm")
pdf(file = paste0("Figures/",species_f,out_base,"_vs_nonspatial.pdf"),
    height = 11,
    width = 8.5)
print(tmap)
dev.off()








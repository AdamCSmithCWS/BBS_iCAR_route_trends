## building a BYM route-level trend model for the BBS
setwd("C:/GitHub/BBS_iCAR_route_trends")
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


firstYear = 2004
lastYear = 2019

scope = "RangeWide"


species = "Blue-headed Vireo"
species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)


out_base <- paste0(species_f,"_Non_spatial_",firstYear,"_",lastYear)


# SPECIES LOOP ------------------------------------------------------------

output_dir <- "output"

csv_files <- dir(output_dir,pattern = out_base,full.names = TRUE)

sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")

    load(sp_file) 
output_dir <- "output"
out_base <- paste0(species_f,"_Non_spatial_",firstYear,"_",lastYear)



stan_data[["N_edges"]] = NULL
stan_data[["node1"]] = NULL
stan_data[["node2"]] = NULL


mod.file = "models/slope_NONiCAR_route.stan"


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
                             sdbeta_rand = runif(1,0.01,0.1),
                             beta_raw_rand = rnorm(stan_data$nroutes,0,0.01))}


slope_stanfit <- slope_model$sample(
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


# export to csv and read in as rstan --------------------------------------













# 


# PLOTTING and trend output -----------------------------------------------

# library(tidybayes)


# route_trajectories <- FALSE #set to FALSE to speed up mapping
# 
# maps = vector(mode = "list",length = 400)
# maps2 = vector(mode = "list",length = 400)
# 
# maps3 = vector(mode = "list",length = 400)
# 
# maps_rand = vector(mode = "list",length = 400)
# maps_space = vector(mode = "list",length = 400)
# 
# trends_out_space <- NULL
# trends_out_rand <- NULL
# sdbeta_dif <- NULL
# sdbeta_space_rand <- NULL

jj <- 0

LC = 0.05
UC = 0.95

output_dir <- "G:/BBS_iCAR_route_trends/output"

strata_map  <- get_basemap(strata_type = strat,
                           transform_laea = TRUE,
                           append_area_weights = FALSE)



loo_stan_spatial <- vector(mode = "list",length = length(selSpecies))
loo_stan_NONspatial <- loo_stan_spatial
trends_out <- NULL


# fix this section to compare spatial and non-spatial ---------------------

library(loo)

for(species in selSpecies){
  
  
  species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)
  
  sp_file_new <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_NONiCAR.RData")
  
  sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")
  
  
  #sp_file <- paste0("output/",species,"Canadian_",firstYear,"_",lastYear,"_slope_route_iCAR2.RData")
  if(file.exists(sp_file_new)){

    load(sp_file)

    sl_rstan <- rstan::read_stan_csv(csv_files)
    
    loo_stan_spatial[[species]] = loo(sl_rstan)
    
    
    # add trend and abundance ----------------------------------------
    
    beta_samples = posterior_samples(sl_rstan,"beta",
                                     dims = "s")
    slopes = beta_samples %>% group_by(s) %>% 
      summarise(b_spat = mean(.value),
                lci_spat = quantile(.value,LC),
                uci_spat = quantile(.value,UC),
                sd_spat = sd(.value),
                prec_spat = 1/var(.value),
                trend_spat = mean((exp(.value)-1)*100),
                lci_trend_spat = quantile((exp(.value)-1)*100,LC),
                uci_trend_spat = quantile((exp(.value)-1)*100,UC),
                .groups = "keep")
    
    alpha_samples = posterior_samples(sl_rstan,"alpha",
                                      dims = "s")
    interc = alpha_samples %>% group_by(s) %>% 
      summarise(abund_spat = mean(exp(.value)),
                lci_i_spat = quantile(exp(.value),LC),
                uci_i_spat = quantile(exp(.value),UC),
                sd_i_spat = sd(exp(.value)),
                prec_i_spat = 1/var(exp(.value)),
                .groups = "keep")
    
    #plot(log(interc$i),slopes$b)
    slops_int = inner_join(slopes,interc,by = "s")
    slops_int$routeF = slops_int$s
    
    route_map_df = as.data.frame(route_map)
    
    route_map_out = left_join(route_map_df,slops_int,by = "routeF")
    route_map_out$species <- species
  

    
    

# non-spatial version -----------------------------------------------------

    

    load(sp_file_new)

    sl_rstan <- rstan::read_stan_csv(csv_files)

    loo_stan_NONspatial[[species]] = loo(sl_rstan)

    # add trend and abundance ----------------------------------------
    
    beta_samples = posterior_samples(sl_rstan,"beta",
                                     dims = "s")
    slopes = beta_samples %>% group_by(s) %>% 
      summarise(b_NONspat = mean(.value),
                lci_NONspat = quantile(.value,LC),
                uci_NONspat = quantile(.value,UC),
                sd_NONspat = sd(.value),
                prec_NONspat = 1/var(.value),
                trend_NONspat = mean((exp(.value)-1)*100),
                lci_trend_NONspat = quantile((exp(.value)-1)*100,LC),
                uci_trend_NONspat = quantile((exp(.value)-1)*100,UC),
                .groups = "keep")
    
    alpha_samples = posterior_samples(sl_rstan,"alpha",
                                      dims = "s")
    interc = alpha_samples %>% group_by(s) %>% 
      summarise(abund_NONspat = mean(exp(.value)),
                lci_i_NONspat = quantile(exp(.value),LC),
                uci_i_NONspat = quantile(exp(.value),UC),
                sd_i_NONspat = sd(exp(.value)),
                prec_i_NONspat = 1/var(exp(.value)),
                .groups = "keep")
    
    #plot(log(interc$i),slopes$b)
    slops_int = inner_join(slopes,interc,by = "s")
    slops_int$routeF = slops_int$s
    
    route_map_out = left_join(route_map_out,slops_int,by = "routeF")
    
    
    
    
    trends_out <- bind_rows(trends_out,route_map_out)
    
    
    
    
  }
  
  
}
  

save(list = c("loo_stan_NONspatial","loo_stan_spatial",
              "trends_out"),
     file = "output/loo_prediction_comparison_spatial_NON.RData")
  


load("output/loo_prediction_comparison_spatial_NON.RData")

# full loo ----------------------------------------------------------------
loo_full <- NULL



for(species in selSpecies){
  
  tmp <- data.frame(looic = loo_stan_spatial[[species]]$estimates["looic","Estimate"],
                    se_looic = loo_stan_spatial[[species]]$estimates["looic","SE"],
                    model = "Spatial")
  
  tmp2 <- data.frame(looic = loo_stan_NONspatial[[species]]$estimates["looic","Estimate"],
                    se_looic = loo_stan_NONspatial[[species]]$estimates["looic","SE"],
                    model = "Non_Spatial")
  
  tmp3 = bind_rows(tmp,tmp2)
  tmp3$species = species
  
  loo_full <- bind_rows(loo_full,tmp3)
  
  }

loo_full$lci <- loo_full$looic - (loo_full$se_looic*(1.96*2))
loo_full$uci <- loo_full$looic + (loo_full$se_looic*(1.96*2))

loo_comp_plot <- ggplot(data = loo_full,aes(x = species,y = looic,colour = model,group = model))+
  geom_point(position = position_dodge(width = 0.3))+
  geom_errorbar(aes(ymin = lci,ymax = uci),alpha = 0.3,width = 0,position = position_dodge(width = 0.3))+
  coord_flip()

print(loo_comp_plot)





# pointwise loo -----------------------------------------------------------
loo_point <- NULL



for(species in selSpecies){
  species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)
  
  
  sp_file <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_iCAR.RData")
  load(sp_file)
  orig_dat_spatial <- data.frame(count_spat = stan_data$count,
                                 year_spat = stan_data$year,
                                 observer_spat = stan_data$observer,
                                 routeF_spat = stan_data$route)
  
  sp_route_file <- paste0(species_f,"_route_data.RData")
  load(paste0("route_maps/",sp_route_file))
  route_df <- as.data.frame(route_map) %>% 
    select(route,routeF,strat)
  
  
  orig_dat_spatial <- orig_dat_spatial %>% 
    left_join(.,route_df, by = c("routeF_spat"="routeF"))
  
  
  # sp_file_new <- paste0(output_dir,"/",species_f,"_",scope,"_",firstYear,"_",lastYear,"_slope_route_NONiCAR.RData")
  # 
  # load(sp_file_new)
  # 
  # 
  # orig_dat_non <- data.frame(count_non = stan_data$count,
  #                                year_non = stan_data$year,
  #                                observer_non = stan_data$observer,
  #                                routeF_non = stan_data$route)
  # 
  # 
  # tmp = bind_cols(orig_dat_spatial,orig_dat_non)
  # 
  
  
  tmp <- as.data.frame(loo_stan_spatial[[species]]$pointwise)
  names(tmp) <- paste(names(tmp),"spatial",sep = "_")
  
  
  tmp2 <- as.data.frame(loo_stan_NONspatial[[species]]$pointwise)
  names(tmp2) <- paste(names(tmp2),"NONspatial",sep = "_")
  
  
  orig_dat <- orig_dat_spatial %>% 
    bind_cols(.,tmp) %>% 
    bind_cols(.,tmp2) %>% 
    mutate(species = species,
           dif_elpd_spat_non = elpd_loo_spatial - elpd_loo_NONspatial)

  loo_point <- bind_rows(loo_point,orig_dat)
  
}


route_n <- loo_point %>% 
  group_by(species,route,routeF_spat) %>% 
  summarise(n_surv = n(),
            mean_count = mean(count_spat),
            mean_dif_elpd = mean(dif_elpd_spat_non))

#obs_mean <- 

loo_point <- loo_point %>% 
  left_join(.,route_n,by = c("species","routeF_spat","route")) %>% 
  mutate(log_count1 = log(count_spat+1))

bi_plot <- ggplot(data = loo_point,aes(x = elpd_loo_spatial,y = elpd_loo_NONspatial))+
  geom_point(aes(colour = n_surv),alpha = 0.2)+
  geom_abline(intercept = 0,slope = 1)+
  facet_wrap(~species,nrow = 3,ncol = 3,scales = "free")


print(bi_plot)



dif_plot_n <- ggplot(data = loo_point,aes(y = dif_elpd_spat_non,x = n_surv))+
  geom_point(alpha = 0.2,position = position_jitter(width = 0.3))+
  geom_abline(intercept = 0,slope = 0)+
  geom_smooth(method = "lm")+
  facet_wrap(~species,nrow = 3,ncol = 3,scales = "free")


print(dif_plot_n)


dif_plot_count <- ggplot(data = loo_point,aes(y = dif_elpd_spat_non,x = log_count1))+
  geom_point(alpha = 0.2,position = position_jitter(width = 0.3))+
  geom_abline(intercept = 0,slope = 0)+
  geom_smooth(method = "lm")+
  facet_wrap(~species,nrow = 3,ncol = 3,scales = "free")


print(dif_plot_count)



dif_plot_mean_count <- ggplot(data = route_n,aes(y = mean_dif_elpd,x = mean_count))+
  geom_point(alpha = 0.2,position = position_jitter(width = 0.3))+
  geom_abline(intercept = 0,slope = 0)+
  scale_x_log10()+
  geom_smooth(method = "lm")+
  facet_wrap(~species,nrow = 3,ncol = 3,scales = "free")


print(dif_plot_mean_count)

# compare trend predictions -----------------------------------------


#### 

trends_comp <- left_join(trends_out,sp_type,by = "species") %>% 
  arrange(.,spatial_pattern,species,routeF) %>% 
  filter(species %in% c("Ring-billed Gull",
                        "Northern Bobwhite")) %>% 
  mutate(cv_abund_spat = sd_spat/abund_spat,
         cv_abund_NONspat = sd_NONspat/abund_NONspat) %>% 
  group_by(species) %>% 
  mutate(rel_prec_trend = scale(sd_NONspat),
         rel_prec_abund = scale(sd_i_NONspat))

bi_plot_trends <- ggplot(data = trends_comp,aes(x = trend_spat,y = trend_NONspat))+
  geom_point(aes(colour = rel_prec_trend),alpha = 0.5)+
  geom_abline(intercept = 0,slope = 1)+
  xlab("Spatial Trend Estimate (%/year)")+
  ylab("Non Spatial Trend Estimate (%/year)")+
  scale_colour_viridis_c(begin = 0, end = 0.9,direction = 1,
                         guide = guide_legend(title = "relative SD of trend",
                                              reverse = TRUE))+
  theme_minimal()+
  #theme(legend.position = "none")+
  facet_wrap(~species,nrow = 2,ncol = 2,scales = "free")

pdf(file = paste0("Figures/bi_plot_trends_spatial_non_by_pred.pdf"),
    width = 8,
    height = 4)
print(bi_plot_trends)
dev.off()


bi_plot_prec <- ggplot(data = trends_comp,aes(x = prec_spat,y = prec_NONspat))+
  geom_point(alpha = 0.3)+
  #geom_point(aes(colour = strat),alpha = 0.5)+
  geom_abline(intercept = 0,slope = 1)+
  xlab("Spatial Precision Trend Estimate (%/year)")+
  ylab("Non Spatial Precision Trend Estimate (%/year)")+
  # scale_colour_viridis_d(begin = 0, end = 0.9,direction = 1,
  #                        guide = guide_legend(title = "relative SD of trend",
  #                                             reverse = TRUE))+
  theme_minimal()+
  theme(legend.position = "none")+
  facet_wrap(~species,nrow = 2,ncol = 2,scales = "free")

print(bi_plot_prec)



# comparing abundance predictions -----------------------------------------



bi_plot_abund <- ggplot(data = trends_comp,aes(x = abund_spat,y = abund_NONspat))+
  geom_point(aes(colour = rel_prec_abund),alpha = 0.5)+
  geom_abline(intercept = 0,slope = 1)+
  xlab("Spatial Abundance Estimate (%/year)")+
  ylab("Non Spatial Abundance Estimate (%/year)")+
  scale_colour_viridis_c(begin = 0, end = 0.9,direction = 1,
                         guide = guide_legend(title = "relative SD of abundance",
                                              reverse = TRUE))+
  theme_minimal()+
  facet_wrap(~species,nrow = 2,ncol = 2,scales = "free")

print(bi_plot_abund)




bi_plot_prec_abund <- ggplot(data = trends_comp,aes(x = prec_i_spat,y = prec_i_NONspat))+
  geom_point(alpha = 0.3)+
  #geom_point(aes(colour = strat),alpha = 0.5)+
  geom_abline(intercept = 0,slope = 1)+
  xlab("Spatial Precision Trend Estimate (%/year)")+
  ylab("Non Spatial Precision Trend Estimate (%/year)")+
  # scale_colour_viridis_d(begin = 0, end = 0.9,direction = 1,
  #                        guide = guide_legend(title = "relative SD of trend",
  #                                             reverse = TRUE))+
  theme_minimal()+
  facet_wrap(~species,nrow = 2,ncol = 2,scales = "free")

print(bi_plot_prec_abund)



bi_plot_prec_abund <- ggplot(data = trends_comp,aes(x = cv_abund_spat,y = cv_abund_NONspat))+
  geom_point(alpha = 0.3)+
  #geom_point(aes(colour = strat),alpha = 0.5)+
  geom_abline(intercept = 0,slope = 1)+
  xlab("Spatial Precision Trend Estimate (%/year)")+
  ylab("Non Spatial Precision Trend Estimate (%/year)")+
  # scale_colour_viridis_d(begin = 0, end = 0.9,direction = 1,
  #                        guide = guide_legend(title = "relative SD of trend",
  #                                             reverse = TRUE))+
  theme_minimal()+
  facet_wrap(~species,nrow = 2,ncol = 2,scales = "free")

print(bi_plot_prec_abund)







    
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

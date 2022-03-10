setwd("C:/GitHub/BBS_iCAR_route_trends/")
library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(sf)
library(spdep)
library(ggforce)
library(brms)
library(patchwork)
source("functions/neighbours_define.R") ## function to define neighbourhood relationships
source("functions/prepare-jags-data-alt.R") ## small alteration of the bbsBayes function
source("functions/get_basemap_function.R") ## loads one of the bbsBayes strata maps
source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output
source("functions/initial_value_functions.R")
## changes captured in a commit on Nov 20, 2020


# load and stratify CASW data ---------------------------------------------
#species = "Pacific Wren"
#species = "Barred Owl"
strat = "bbs_usgs"
model = "slope"
scope = "RangeWide"

firstYear = 2004
lastYear = 2019 # final year to consider

# select a minimum year for prediction (i.e., a route has to have data between 2004 and 2011 to be included)
# similar to "L" in Burkner et al 2020 (https://doi.org/10.1080/00949655.2020.1783262)
minimumYear = 2011 

species = "Grasshopper Sparrow"
species = "Blue Jay"
species = "Blue-headed Vireo"
species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)



cv_sum = NULL

for(sppn in c("iCAR","BYM","Non_spatial")){
  
  load(paste0("data/",species_f,"CV_base_data.RData"))
  
  output_dir <- "output"
  spp <- paste0("_",sppn,"_")
 load(paste0("output/",species_f,spp,"_pred_save.RData")) 
 
 cv_sum = bind_rows(cv_sum,predictions_save)
 
 
 
 
}

n_obs <- full_obs_df %>% 
  group_by(route,ObsN) %>% 
  summarise(n_years = n()) %>% 
  group_by(route) %>% 
  summarise(sum_n_years = sum(n_years),
            max_nyears = max(n_years),
            mean_nyears = mean(n_years),
            n_obs = n())


# point wise differences among models -------------------------------------
cv_sum$Year <- factor(cv_sum$r_year)

diffs <- cv_sum %>% 
  select(count,Year,r_year,route,observer,E_pred_i,model,log_lik_mean) %>% 
  pivot_wider(.,names_from = model,
              values_from = log_lik_mean) %>% 
  mutate(iCAR_BYM = iCAR - BYM,
         iCAR_Non_spatial = iCAR-Non_spatial,
         BYM_Non_spatial = BYM-Non_spatial) %>% 
  left_join(.,n_obs,by = "route")

cv_sum <- cv_sum %>% 
  left_join(.,n_obs,by = "route") %>% 
  mutate(model = factor(model,levels = c("Non_spatial","BYM","iCAR"),
                        ordered = TRUE))

save(list = c("diffs","cv_sum"),
     file = "data/example_cv_comparisons.RData")
# save(list = c("diffs","cv_sum"),
#      file = "data/GRSP_cv_comparisons.RData")



diffs_p1 <- ggplot(data = diffs,aes(x = count+1,y = iCAR_BYM))+
  geom_point(alpha = 0.2)+
  scale_x_continuous(trans = "log10")+
  geom_smooth()
diffs_p2 <- ggplot(data = diffs,aes(x = count+1,y = iCAR_Non_spatial))+
  geom_point(alpha = 0.2)+
  scale_x_continuous(trans = "log10")+
  geom_smooth()
diffs_p3 <- ggplot(data = diffs,aes(x = count+1,y = BYM_Non_spatial))+
  geom_point(alpha = 0.2)+
  scale_x_continuous(trans = "log10")+
  geom_smooth()
print(diffs_p1 + diffs_p2 + diffs_p3)


diffs_p1 <- ggplot(data = diffs,aes(x = sum_n_years,y = iCAR_BYM))+
  geom_point(alpha = 0.2)+
  geom_smooth()
diffs_p2 <- ggplot(data = diffs,aes(x = sum_n_years,y = iCAR_Non_spatial))+
  geom_point(alpha = 0.2)+
  geom_smooth()
diffs_p3 <- ggplot(data = diffs,aes(x = sum_n_years,y = BYM_Non_spatial))+
  geom_point(alpha = 0.2)+
  geom_smooth()
print(diffs_p1 + diffs_p2 + diffs_p3)


diffs_p1 <- ggplot(data = diffs,aes(x = max_nyears,y = iCAR_BYM))+
  geom_point(alpha = 0.2)+
  geom_smooth()
diffs_p2 <- ggplot(data = diffs,aes(x = max_nyears,y = iCAR_Non_spatial))+
  geom_point(alpha = 0.2)+
  geom_smooth()
diffs_p3 <- ggplot(data = diffs,aes(x = max_nyears,y = BYM_Non_spatial))+
  geom_point(alpha = 0.2)+
  geom_smooth()
print(diffs_p1 + diffs_p2 + diffs_p3)






lpos = function(x){
  p = length(which(x > 0))/length(x)
}
mndiffs = diffs %>% 
  group_by(Year) %>% 
  summarise(m_iCAR_BYM = median(iCAR_BYM),
            m_iCAR_Non_spatial = median(iCAR_Non_spatial),
            m_BYM_Non_spatial = median(BYM_Non_spatial),
            nbet_iCAR_BYM = lpos(iCAR_BYM),
            nbet_iCAR_Non_spatial = lpos(iCAR_Non_spatial),
            nbet_BYM_Non_spatial = lpos(BYM_Non_spatial))
mndiffs
mndiffs = diffs %>% 
  summarise(m_iCAR_BYM = median(iCAR_BYM),
            m_iCAR_Non_spatial = median(iCAR_Non_spatial),
            m_BYM_Non_spatial = median(BYM_Non_spatial),
            nbet_iCAR_BYM = lpos(iCAR_BYM),
            nbet_iCAR_Non_spatial = lpos(iCAR_Non_spatial),
            nbet_BYM_Non_spatial = lpos(BYM_Non_spatial))
mndiffs

spat_non = ggplot(data = diffs,aes(y = iCAR_Non_spatial,x = BYM_Non_spatial))+
  geom_point(alpha = 0.05)+
  geom_abline(slope = 1,intercept = 0)
spat_nonzoom = ggplot(data = diffs,aes(y = iCAR_Non_spatial,x = BYM_Non_spatial))+
  geom_point(alpha = 0.05)+
  geom_abline(slope = 1,intercept = 0)+
  coord_cartesian(xlim = c(-2,2),
                  ylim = c(-2,2))

print(spat_non + spat_nonzoom)


m1 = brm(log_lik_mean ~ model*max_nyears + (1|E_pred_i),
         data = cv_sum,
         save_pars = save_pars(all = TRUE))


m2 = brm(log_lik_mean ~ model + max_nyears + (1|E_pred_i),
         data = cv_sum,
         save_pars = save_pars(all = TRUE))


summary(m1)
summary(m2)
lm1 = loo(m1,moment_match = TRUE)
lm2 = loo(m2,moment_match = TRUE)



  
setwd("C:/GitHub/BBS_iCAR_route_trends/")
library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(sf)
library(spdep)
library(ggforce)
library(brms)
library(patchwork)
library(lme4)
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



load("Data/sp_sel.RData")

cv_sum = NULL
n_obs <- NULL
for(species in sp_sel[1:5]){
  
  species_f <- gsub(species,pattern = " ",replacement = "_",fixed = T)
  species_f <- gsub(species_f,pattern = "'",replacement = "",fixed = T)
  



for(sppn in c("iCAR","BYM","Non_spatial")){
  
  load(paste0("data/",species_f,"CV_base_data.RData"))
  
  output_dir <- "output"
  spp <- paste0("_",sppn,"_")
  load(paste0("output/",species_f,spp,"_pred_save.RData")) 

  cv_sum = bind_rows(cv_sum,predictions_save)
  
  
  
  
}
  n_obst <- full_obs_df %>% 
    group_by(route,ObsN) %>% 
    summarise(n_years = n()) %>% 
    group_by(route) %>% 
    summarise(sum_n_years = sum(n_years),
              max_nyears = max(n_years),
              mean_nyears = mean(n_years),
              n_obs = n(),
              .groups = "keep") %>% 
    mutate(species = species)
  
  n_obs <- bind_rows(n_obs,n_obst)
  
}




# point wise differences among models -------------------------------------
cv_sum$Year <- factor(cv_sum$r_year)

diffs <- cv_sum %>% 
  select(species,count,Year,r_year,route,observer,E_pred_i,model,log_lik_mean) %>% 
  pivot_wider(.,names_from = model,
              values_from = log_lik_mean) %>% 
  mutate(iCAR_BYM = iCAR - BYM,
         iCAR_Non_spatial = iCAR-Non_spatial,
         BYM_Non_spatial = BYM-Non_spatial) %>% 
  left_join(.,n_obs,by = c("route","species"))

cv_sum <- cv_sum %>% 
  left_join(.,n_obs,by = c("route","species")) %>% 
  mutate(model = factor(model,levels = c("Non_spatial","BYM","iCAR"),
                        ordered = FALSE),
         I = paste(species,E_pred_i,sep = "-"))

save(list = c("diffs","cv_sum"),
     file = "data/cv_summary_25_data.RData")
# save(list = c("diffs","cv_sum"),
#      file = "data/GRSP_cv_comparisons.RData")


# 
# diffs_p1 <- ggplot(data = diffs,aes(x = count+1,y = iCAR_BYM))+
#   geom_point(alpha = 0.2)+
#   scale_x_continuous(trans = "log10")+
#   geom_smooth()
# diffs_p2 <- ggplot(data = diffs,aes(x = count+1,y = iCAR_Non_spatial))+
#   geom_point(alpha = 0.2)+
#   scale_x_continuous(trans = "log10")+
#   geom_smooth()
# diffs_p3 <- ggplot(data = diffs,aes(x = count+1,y = BYM_Non_spatial))+
#   geom_point(alpha = 0.2)+
#   scale_x_continuous(trans = "log10")+
#   geom_smooth()
# print(diffs_p1 + diffs_p2 + diffs_p3)
# 
# 
# diffs_p1 <- ggplot(data = diffs,aes(x = sum_n_years,y = iCAR_BYM))+
#   geom_point(alpha = 0.2)+
#   geom_smooth()
# diffs_p2 <- ggplot(data = diffs,aes(x = sum_n_years,y = iCAR_Non_spatial))+
#   geom_point(alpha = 0.2)+
#   geom_smooth()
# diffs_p3 <- ggplot(data = diffs,aes(x = sum_n_years,y = BYM_Non_spatial))+
#   geom_point(alpha = 0.2)+
#   geom_smooth()
# print(diffs_p1 + diffs_p2 + diffs_p3)
# 
# 
# diffs_p1 <- ggplot(data = diffs,aes(x = max_nyears,y = iCAR_BYM))+
#   geom_point(alpha = 0.2)+
#   geom_smooth()
# diffs_p2 <- ggplot(data = diffs,aes(x = max_nyears,y = iCAR_Non_spatial))+
#   geom_point(alpha = 0.2)+
#   geom_smooth()
# diffs_p3 <- ggplot(data = diffs,aes(x = max_nyears,y = BYM_Non_spatial))+
#   geom_point(alpha = 0.2)+
#   geom_smooth()
# print(diffs_p1 + diffs_p2 + diffs_p3)
# 





lpos = function(x){
  p = length(which(x > 0))/length(x)
}
mndiffs = diffs %>% 
  group_by(Year,species) %>% 
  summarise(m_iCAR_BYM = median(iCAR_BYM),
            m_iCAR_Non_spatial = median(iCAR_Non_spatial),
            m_BYM_Non_spatial = median(BYM_Non_spatial),
            nbet_iCAR_BYM = lpos(iCAR_BYM),
            nbet_iCAR_Non_spatial = lpos(iCAR_Non_spatial),
            nbet_BYM_Non_spatial = lpos(BYM_Non_spatial))
mndiffs
mndiffs = diffs %>% 
  group_by(species) %>% 
  summarise(m_iCAR_BYM = median(iCAR_BYM),
            m_iCAR_Non_spatial = median(iCAR_Non_spatial),
            m_BYM_Non_spatial = median(BYM_Non_spatial),
            nbet_iCAR_BYM = lpos(iCAR_BYM),
            nbet_iCAR_Non_spatial = lpos(iCAR_Non_spatial),
            nbet_BYM_Non_spatial = lpos(BYM_Non_spatial))
mndiffs

# spat_non = ggplot(data = diffs,aes(y = iCAR_Non_spatial,x = BYM_Non_spatial))+
#   geom_point(alpha = 0.05)+
#   geom_abline(slope = 1,intercept = 0)
# spat_nonzoom = ggplot(data = diffs,aes(y = iCAR_Non_spatial,x = BYM_Non_spatial))+
#   geom_point(alpha = 0.05)+
#   geom_abline(slope = 1,intercept = 0)+
#   coord_cartesian(xlim = c(-2,2),
#                   ylim = c(-2,2))
# 
# print(spat_non + spat_nonzoom)

# 
# m1 = brm(log_lik_mean ~ model*max_nyears + (1|E_pred_i),
#          data = cv_sum,
#          save_pars = save_pars(all = TRUE))
# 
# 
 # m2 = brm(log_lik_mean ~ model*species + max_nyears + (1|I),
 #          data = cv_sum,
 #          save_pars = save_pars(all = TRUE))
 
 # m2a = lmer(log_lik_mean ~ model + species + max_nyears + (1|I),
 #            data = cv_sum)
 
 m3a = lm(iCAR_Non_spatial ~ species + max_nyears,
            data = diffs)
 
 m4a = lm(BYM_Non_spatial ~ species + max_nyears,
          data = diffs)
 m5a = lm(iCAR_BYM ~ species + max_nyears,
          data = diffs)
 
 
summary(m3a)

newdat = expand.grid(species = unique(diffs$species),
                     max_nyears = 10)

predictn4 <- predict(m4a,newdata = newdat,se.fit = TRUE) %>% 
  as.data.frame() %>% 
  mutate(comparison = "BYM_Non_spatial",
         lci = fit-(se.fit*1.96),
         uci = fit+(se.fit*1.96),
         species = newdat$species)
predictn3 <- predict(m3a,newdata = newdat,se.fit = TRUE) %>% 
  as.data.frame() %>% 
  mutate(comparison = "iCAR_Non_spatial",
         lci = fit-(se.fit*1.96),
         uci = fit+(se.fit*1.96),
         species = newdat$species)

predictn5 <- predict(m5a,newdata = newdat,se.fit = TRUE) %>% 
  as.data.frame() %>% 
  mutate(comparison = "iCAR_BYM",
         lci = fit-(se.fit*1.96),
         uci = fit+(se.fit*1.96),
         species = newdat$species)

predicts <- bind_rows(predictn3,predictn4,predictn5)

fig1 <- ggplot(data = predicts,aes(x = species,y = fit,
                                   group = comparison,
                                   colour = comparison),
               position = position_dodge(width = 0.2))+
  geom_errorbar(aes(ymin = lci,ymax = uci),
                alpha = 0.3,
                width = 0,
                position = position_dodge(width = 0.2))+
  geom_point(position = position_dodge(width = 0.2))+
  theme_bw()+
  labs(subtitle = "Positive values support the more spatial model")+
  xlab("")+
  ylab("Mean pointwise difference in posterior log-likelihood")+
  coord_flip(ylim = c(-1,1))

pdf(file = "Figures/temporary_comparison_summary.pdf")
print(fig1)
dev.off()



# 
# 
# summary(m1)
 summary(m2a)

 newdat = expand.grid(model = levels(cv_sum$model),
                     species = unique(cv_sum$species),
                     max_nyears = 10)

 predictn <- predict(m2a,newdata = newdat,re.form = NA) %>% 
   as.data.frame()
   
 newdat <- bind_cols(newdat,predictn)

 
 fig1 <- ggplot(data = predictn,aes(x = species,y = Estimate,
                                    group = model, colour = model))+
   geom_errorbar(aes(ymin = lci,ymax = uci),
                 alpha = 0.3)+
   geom_point()+
   theme_bw()+
   coord_flip(xlim = c(-1,1))
 
 print(fig1)
 
 
 
 


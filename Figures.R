# Figures to support the paper
library(ggplot2)
library(sf)


# 1 spatial set-up demo -----------------------------------------------------

load(paste0("Data/Blue-headed_Vireo_RangeWide_route_maps_data.RData"))

fig1 <- ggplot(data = centres)+ 
  geom_sf(aes(col = strat_lab)) + 
  geom_segment(data=DA,aes(x = long, y = lat,xend=long_to,yend=lat_to),inherit.aes = FALSE,size=0.3,alpha=0.1) +
  geom_sf(data = vintj,alpha = 0,colour = grey(0.95))+ 
  geom_sf(data = real_strata_map,alpha = 0,colour = grey(0.85))

print(fig1)

f = 1
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)

# 2 trend maps comparing 3 models -----------------------------------------


f = 2
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)

# 3 Trend correlations among models ------------------------------------------------------


f = 3
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)

# 4 ABundance correlations among models -----------------------------------


f = 4
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)

# ?5 Trajectories??? ------------------------------------------------------



f = "5alt"
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)

# 5 Cross validation results for example species --------------------------



f = 5
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)

# 6 Cross validation results for ~25 species ------------------------------



f = 6
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)

# 7 trend maps for other species ------------------------------------------




f = 7
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)

# 8 Group level mean trend maps -------------------------------------------




f = 7
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)










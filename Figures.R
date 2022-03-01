# Figures to support the paper
library(spdep)
library(sf)
library(tidyverse)


# 1 spatial set-up demo -----------------------------------------------------

load(paste0("Data/Blue-headed_Vireo_RangeWide_route_maps_data.RData"))

base_map = bbsBayes::load_map(stratify_by = "bbs_usgs")
box <- st_coordinates(st_as_sfc(st_bbox(centres)))
xl <- range(box[,"X"])
yl <- range(box[,"Y"])

fig <- ggplot(data = centres)+ 
  geom_sf(data = base_map,alpha = 0,colour = grey(0.9),size = 0.2)+
  geom_sf(size = 0.5) + 
  geom_segment(data=DA,aes(x = long, y = lat,xend=long_to,yend=lat_to),
               inherit.aes = FALSE,size=0.3,alpha=0.1) +
  #geom_sf(data = vintj,alpha = 0,colour = grey(0.85))+
  theme_minimal() +
  xlab("")+
  ylab("")+
  coord_sf(xlim = xl,
           ylim = yl)+
  theme(legend.position = "none")

#print(fig)

f = 1
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)
print(fig)
dev.off()

# 2 trend maps comparing 3 models -----------------------------------------



f = 2
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)
print(fig)
dev.off()

# 3 Trend correlations among models ------------------------------------------------------


f = 3
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)
print(fig)
dev.off()

# 4 ABundance correlations among models -----------------------------------


f = 4
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)
print(fig)
dev.off()

# ?5 Trajectories??? ------------------------------------------------------



f = "5alt"
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)
print(fig)
dev.off()

# 5 Cross validation results for example species --------------------------



f = 5
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)
print(fig)
dev.off()

# 6 Cross validation results for ~25 species ------------------------------



f = 6
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)
print(fig)
dev.off()

# 7 trend maps for other species ------------------------------------------




f = 7
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)
print(fig)
dev.off()

# 8 Group level mean trend maps -------------------------------------------




f = 7
h = 5
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)
print(fig)
dev.off()










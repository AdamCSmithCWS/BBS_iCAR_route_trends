# Figures to support the paper
library(spdep)
library(sf)
library(tidyverse)
library(patchwork)


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


load("Figures/example_trend_map_data.RData")

strata_map  <- bbsBayes::load_map(stratify_by = "bbs_usgs")


strata_bounds <- st_union(trend_plot_map) #union to provide a simple border of the realised strata
bb = st_bbox(strata_bounds)
xlms = as.numeric(c(bb$xmin,bb$xmax))
ylms = as.numeric(c(bb$ymin,bb$ymax))


breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
lgnd_head <- "Trend\n"
trend_title <- "Trend"
labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
labls = paste0(labls, " %/year")
trend_plot_map$Tplot <- cut(trend_plot_map$trend,breaks = c(-Inf, breaks, Inf),labels = labls)

map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                 "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
names(map_palette) <- labls


fig = ggplot(trend_plot_map)+
  #geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
  geom_sf(aes(colour = Tplot,size = abund))+
  scale_size_continuous(range = c(0.05,2),
                        name = "Mean Abundance")+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0(lgnd_head,2004,"-",2019))+
  coord_sf(xlim = xlms,ylim = ylms)+
  facet_wrap(vars(version),nrow = 3)+
  theme_bw()


f = 2
h = 10
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)
print(fig)
dev.off()

# 3 Trend correlations among models ------------------------------------------------------
# iCAR vs non-spatial trends and iCAR vs BYM trends

load("Figures/example_trend_comparison_data.RData")


tcplot <- ggplot(data = trend_comp,
                 aes(x = trend_BYM,y = trend_Non_spatial,colour = trend_se_BYM))+
  geom_abline(slope = 1,intercept = 0)+
  geom_point(alpha = 0.5,size = 0.75)+
  ylab("Trend non-spatial model")+
  xlab("Trend BYM model")+
  scale_colour_viridis_c(aesthetics = "colour",direction = 1)+
  theme_classic()+
  scale_y_continuous(limits = c(-13,13))+
  scale_x_continuous(limits = c(-13,13))+
  theme(legend.position = "none")

tcplot2 <- ggplot(data = trend_comp,
                  aes(y = trend_iCAR,x = trend_BYM,colour = trend_se_BYM))+
  geom_abline(slope = 1,intercept = 0)+
  geom_point(alpha = 0.5,size = 0.75)+
  ylab("Trend iCAR model")+
  xlab("Trend BYM model")+
  scale_colour_viridis_c(aesthetics = "colour",direction = 1)+
  theme_classic()+
  scale_y_continuous(limits = c(-13,13))+
  scale_x_continuous(limits = c(-13,13))+
  theme(legend.position = "none")


fig <- (tcplot + tcplot2)+
  plot_layout(ncol = 2,
              guides = "collect")



f = 3
h = 4
w = 7

pdf(file = paste0("Figures/Figure_",f,".pdf"),
    height = h,
    width = w)
print(fig)
dev.off()

# 4 ABundance correlations among models -----------------------------------


tcplot <- ggplot(data = abund_comp,
                 aes(x = abund_BYM,y = abund_Non_spatial,colour = abund_se_BYM))+
  geom_abline(slope = 1,intercept = 0)+
  geom_point(alpha = 0.5,size = 0.75)+
  ylab("abund non-spatial model")+
  xlab("abund BYM model")+
  scale_colour_viridis_c(aesthetics = "colour",direction = 1)+
  theme_classic()+
  scale_y_continuous(limits = c(0,13))+
  scale_x_continuous(limits = c(0,13))+
  theme(legend.position = "none")

tcplot2 <- ggplot(data = abund_comp,
                  aes(y = abund_iCAR,x = abund_BYM,colour = abund_se_BYM))+
  geom_abline(slope = 1,intercept = 0)+
  geom_point(alpha = 0.5,size = 0.75)+
  ylab("abund iCAR model")+
  xlab("abund BYM model")+
  scale_colour_viridis_c(aesthetics = "colour",direction = 1)+
  theme_classic()+
  scale_y_continuous(limits = c(0,13))+
  scale_x_continuous(limits = c(0,13))+
  theme(legend.position = "none")


fig <- (tcplot + tcplot2)+
  plot_layout(ncol = 2,
              guides = "collect")

f = 4
h = 4
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


load("data/example_cv_comparisons.RData")

#re-orient, adjust x-labels, include count of points on either side

## or replace with trend-biplot of two spatial models
diffs_p1 <- ggplot(data = diffs,aes(x = max_nyears,y = iCAR_BYM))+
  geom_point(alpha = 0.2)+
  geom_smooth()+
  theme_classic()
diffs_p2 <- ggplot(data = diffs,aes(x = max_nyears,y = iCAR_Non_spatial))+
  geom_point(alpha = 0.2)+
  geom_smooth()+
  theme_classic()+
  xlab("Maximum year span of single observer")
diffs_p3 <- ggplot(data = diffs,aes(x = max_nyears,y = BYM_Non_spatial))+
  geom_point(alpha = 0.2)+
  geom_smooth()+
  theme_classic()
fig <- (diffs_p1 + diffs_p2 + diffs_p3)


comp = ggplot(data = cv_sum,)

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










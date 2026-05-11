library(sf)
library(terra)
library(tidyverse)
library(ggridges)
library(viridis)
# library(hrbrthemes)

library(viridisLite)
library(colorspace)

# Function to make a darker, matte version of viridis-H
viridis_matte <- function(n = 256,
                          option = "H",
                          desat_amount = 0.4,   # 0 = no change, 1 = fully gray
                          darken_amount = 0.2) {# 0 = no change, 1 = black
  # 1) base viridis palette
  pal <- viridisLite::viridis(n, option = option, direction = -1)
  
  # 2) reduce saturation (chroma)
  pal <- colorspace::desaturate(pal, amount = desat_amount)
  
  # 3) darken (reduce luminance/value)
  pal <- colorspace::darken(pal, amount = darken_amount)
  
  pal
}

# World basemap via rnaturalearth (land polygons)
wrld_org <- rnaturalearth::ne_states(country = 'united states of america', 
                                     returnclass = "sf")

ca <- wrld_org %>% filter(name_en == 'California')

#base layer to extend when plants are not in whole California
r <-  rast("figures/osli_all_means.tif") %>% 
  crop(ca, mask=T)
wrld <- wrld_org %>% 
  st_crop(extend(ext(r), 0.5))

wrldp <- st_transform(wrld, crs=4269)

curr <- rast('outputs/models/CA_osli_current_ensemble.tif') %>% mean()/1000

# hist(values(curr))
# hist(values(r[[9]]))
currM <- mask(curr, r[[1]])
dats <- 
  tibble(
    current = values(currM, na.rm=T)[,1] %>% unname(), 
    # ssp245_2044 = values(r[[1]], na.rm=T)[,1] %>% unname(),
    # ssp245_2074 = values(r[[2]], na.rm=T)[,1] %>% unname(),
    ssp245_2100 = values(r[[3]], na.rm=T)[,1] %>% unname(),
    
    # ssp585_2044 = values(r[[7]], na.rm=T)[,1] %>% unname(),
    # ssp585_2074 = values(r[[8]], na.rm=T)[,1] %>% unname(),
    ssp585_2100 = values(r[[9]], na.rm=T)[,1] %>% unname()
  ) %>% pivot_longer(
      cols = everything()
  ) %>% 
  mutate(Scenario = factor(name, 
                           levels= c(
                             # 'ssp245_2044', 
                             # 'ssp245_2074', 
                             'ssp585_2100',
                             'ssp245_2100',
                             'current'
                             # 'ssp585_2044', 
                             # 'ssp585_2074', 
                                     )))

ggplot(dats, aes(x = value, y = Scenario, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 1.2, linewidth = 0.1, rel_min_height = 0.005) +
  scale_fill_gradientn(
    colours = c('#d7191c','#fdae61', '#ffffbf', '#abdda4', '#2b83ba')
    # colours = viridis_matte(256, desat_amount = 0.15, darken_amount = .1),
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by=.15)
    )+
  # scale_fill_viridis(name = "Temp. [F]", option = "H",alpha = 0.85 , direction = -1) +
  # scale_y_continuous(expand = expansion(0.1)) +
  geom_vline(xintercept = 0.665, linetype=2) +
  labs(x = 'Suitability') +
  theme_light(base_size = 18) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(
    legend.position="none",
    # panel.spacing = unit(10, "lines"),
    # strip.text.x = element_text(size = 12)
  )

ggsave('figures/density_suit_draft.png', dpi=600)  
ggsave('figures/density_suit_draft.svg', dpi=600)  



#### PLANTS ####
#base layer to extend when plants are not in whole California
r <-  rast("figures/plants_all_future_rich_raw.tif") %>% 
  crop(ca, mask=T)
wrld <- wrld_org %>% 
  st_crop(extend(ext(r), 0.5))

wrldp <- st_transform(wrld, crs=4269)

curr <- rast('figures/plants_current_richness.tif') %>% resample(r[[1]])

currM <- mask(curr, r[[1]])
names(r)
dats <- 
  tibble(
    current = values(currM, na.rm=T)[,1] %>% unname(), 
    # ssp245_2044 = values(r[[1]], na.rm=T)[,1] %>% unname(),
    # ssp245_2074 = values(r[[2]], na.rm=T)[,1] %>% unname(),
    ssp245_2100 = values(r[[7]], na.rm=T)[,1] %>% unname(),
    
    # ssp585_2044 = values(r[[7]], na.rm=T)[,1] %>% unname(),
    # ssp585_2074 = values(r[[8]], na.rm=T)[,1] %>% unname(),
    ssp585_2100 = values(r[[9]], na.rm=T)[,1] %>% unname()
  ) %>% pivot_longer(
    cols = everything()
  ) %>% 
  mutate(Scenario = factor(name, 
                           levels= c(
                             # 'ssp245_2044', 
                             # 'ssp245_2074', 
                             'ssp585_2100',
                             'ssp245_2100',
                             'current'
                             # 'ssp585_2044', 
                             # 'ssp585_2074', 
                           )))

ggplot(dats, aes(x = value, y = Scenario, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 1, linewidth = 0.1,
                               rel_min_height = 0.005) +
  scale_fill_gradientn(
    colours = viridis_matte(256,option = 'D', desat_amount = 0.15, darken_amount = .1),
  ) +
  # ylim(NA, 'current')+
  labs(x = 'Richness') +
  theme_light(base_size = 14) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position="none")

ggsave('figures/plant_density_suit_draft.png', dpi=600)  
ggsave('figures/plant_density_suit_draft.svg', dpi=600)  

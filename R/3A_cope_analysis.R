library(sf)
library(terra)
library(tidyterra)
library(tidyverse)
library(rnaturalearth)

ca <- rnaturalearth::ne_states(country = 'united states of america') %>% 
  filter(name == 'California')

cope_index <- function(sf, rf, rc, k, v = 0.5, delta = 1, h = 3){
  

  g_rc <- rc / (rc + h)
  
  wr <- 1 + v * g_rc * tanh(k * log((rf + delta)/(rc + delta)))
  
  cope <- sf * wr
  return(cope)
}

sf <- rast('outputs/models/projections/means/ssp585_2075-2100_EMmean.tif') %>% 
  crop(ca, mask=T)/1000

rc <- rast('outputs/baseline_floral_richness.tif')

rf <- rast('outputs/2100_ssp585_floral_floral_richness.tif')

COPE_osmia <- cope_index(sf, rf, rc, k=0.5, delta = 1, h = 3)

# number of species changes
ggplot() + 
  geom_spatraster(data = COPE_osmia) + 
  scale_fill_terrain_c()+
  theme_minimal()

#To plot in QGIS
writeRaster(COPE_osmia %>% 
              disagg(fact=8, method='bilinear') %>%
              mask(ca), 'outputs/COPE_index_osmia.tif', overwrite=T)

gainloss_COPE <- COPE_osmia - sf

ggplot() + 
  geom_spatraster(data = gainloss_COPE) + 
  scale_fill_gradient2(low='darkred', mid='gray90', high='navyblue', na.value = 'transparent')+
  theme_minimal()

#To plot in QGIS
writeRaster(gainloss_COPE %>% 
              disagg(fact=8, method='bilinear') %>%
              mask(ca), 'outputs/gainloss_COPE_index_osmia.tif', overwrite=T)
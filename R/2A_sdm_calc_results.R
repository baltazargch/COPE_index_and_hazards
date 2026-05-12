# libraries
library(sf)
library(terra)
library(tidyterra)
library(tidyverse)
library(patchwork)
library(rnaturalearth)


# Osmia lignaria model results --------------------------------------------

# General Osmia lignaria models performance results
occs_osmia <- read_csv('inputs/records/pres_abs_oslig_target_group.csv')

# n unique records
occs_osmia %>% 
  filter(class != 'background') %>% 
  count()

# results by pseudo-absence number
ens_res_size <- read_csv("outputs/csv/perf_by_PA_size_auc_gated.csv")

# PA1 (2000) PA2 (4000) PA3 (7000) P4 (10000)
ens_res_size %>% arrange(PA)

# final single models metrics
ens_res_selected <- read_csv("outputs/csv/tuning_summary_by_size.csv")

ens_res_selected %>% 
  filter(PA == min(PA, na.rm = TRUE)) %>% 
  pivot_longer(
    cols = TSS_median:KAPPA_median,
    names_to = "metric",
    values_to = "value"
  ) %>% 
  group_by(metric) %>% 
  summarise(
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

# final ensemble model metrics
mods_ens <- readRDS('outputs/models/osmia_bm_ensemble_model.rds')

biomod2::get_evaluations(mods_ens) %>%
  group_by(metric.eval) %>% 
  summarise(
    mean_value = mean(calibration, na.rm = TRUE)
  )

# General plot for exploring model's maps baseline vs future
ca <- rnaturalearth::ne_states(country = 'united states of america') %>% 
  filter(name == 'California')

base_osmia <- rast('outputs/models/NA_osli_current_ensemble.tif') %>% 
  mean()/1000

NA_osmia <- ggplot() + 
  geom_spatraster(data = base_osmia) + 
  scale_fill_viridis_c(option='H', na.value = NA, alpha = 0.8) + 
  theme_minimal()
  
CA_osmia <- ggplot() + 
  geom_spatraster(data = crop(base_osmia, ca, mask=T)) + 
  scale_fill_viridis_c(option='H', na.value = NA, alpha = 0.8) + 
  theme_minimal()  

# descriptive spatial patterns in results. 
NA_osmia / CA_osmia + plot_layout(guides = 'collect')

# raster to plot in QGIS, smooth Fig 1A  
rout <- crop(base_osmia, ca, mask=T) %>% 
  disagg(fact=8, method='bilinear') %>%
  mask(ca)

writeRaster(rout, 'outputs/ca_baseline_osmia.tif')

# prepare Fig 1B - high-emissions only

# Table 1. suitable areas for Osmia lignaria across scenarios


# Floral resources model results ------------------------------------------

# General stats for floral resources models

# Baseline and future richness patterns



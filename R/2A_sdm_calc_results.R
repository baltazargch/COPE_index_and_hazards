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

writeRaster(rout, 'outputs/ca_baseline_osmia.tif', overwrite=TRUE)

# prepare Fig 1B - high-emissions only

ssp585_osmia <- rast('outputs/models/projections/means/ssp585_2075-2100_EMmean.tif')

ca_ssp585_osmia <- crop(ssp585_osmia, ca, mask=T)/1000

pct_change <- (ca_ssp585_osmia - crop(base_osmia, ca, mask=T)) / crop(base_osmia, ca, mask=T) *100

writeRaster(disagg(pct_change, fact=8, method='bilinear') %>% mask(ca), 
            'outputs/pct_change_high_emission_osmia.tif', overwrite=TRUE)

# Table 1. suitable areas for Osmia lignaria across scenarios
# California area in km2
ca_area <- as.numeric(st_area(ca)) * 1e-6

# 10th percentile threshold from training occurrences
cutoff <- terra::extract(
  rout,
  occs_osmia %>% select(lon, lat)
) %>% 
  select(-ID) %>% 
  pull(1) %>% 
  quantile(probs = 0.10, na.rm = TRUE) %>% 
  as.numeric()

# Baseline suitable area
baseline_area <- expanse(
  rout >= cutoff,
  unit = "km",
  byValue = TRUE
) %>% 
  as_tibble() %>% 
  filter(value == 1) %>% 
  pull(area)

# Future rasters
fut_rast <- rast(list.files(
  "outputs/models/projections/means/",
  pattern = "\\.tif$",
  full.names = TRUE
))

names(fut_rast) <- sources(fut_rast) %>% 
  basename() %>% 
  str_remove("_EMmean\\.tif$")

# If future rasters are scaled by 1000
fut_rast <- mask(fut_rast, ca) / 1000

# Function to calculate suitable area for one raster layer
calc_suitable_area <- function(r, cutoff) {
  expanse(r >= cutoff, unit = "km", byValue = TRUE) %>% 
    as_tibble() %>% 
    filter(value == 1) %>% 
    pull(area)
}

# Future table
future_table <- map_dfr(seq_len(nlyr(fut_rast)), \(i) {
  
  layer_name <- names(fut_rast)[i]
  
  tibble(
    layer = layer_name,
    suitable_area_km2 = calc_suitable_area(fut_rast[[i]], cutoff)
  )
})

# Parse scenario and period from layer names
# Modify these regexes if your layer names are different
future_table <- future_table %>% 
  mutate(
    scenario = str_extract(layer, "ssp\\d+|SSP\\d+") %>% str_to_upper(),
    period = case_when(
      str_detect(layer, "2015|2044") ~ "2015--2044",
      str_detect(layer, "2045|2074") ~ "2045--2074",
      str_detect(layer, "2075|2100") ~ "2075--2100",
      TRUE ~ NA_character_
    )
  ) %>% 
  mutate(
    period = factor(
      period,
      levels = c("2015--2044", "2045--2074", "2075--2100")
    )
  )

# Add baseline row and calculate changes
table_osmia_area <- bind_rows(
  tibble(
    layer = "Baseline",
    scenario = "Baseline",
    period = "Baseline",
    suitable_area_km2 = baseline_area
  ),
  future_table
) %>% 
  mutate(
    pct_california = 100 * suitable_area_km2 / ca_area,
    change_from_baseline_km2 = suitable_area_km2 - baseline_area,
    pct_change_from_baseline = 100 * (suitable_area_km2 - baseline_area) / baseline_area
  ) %>% 
  mutate(
    across(
      c(
        suitable_area_km2,
        pct_california,
        change_from_baseline_km2,
        pct_change_from_baseline
      ),
      \(x) round(x, 1)
    )
  ) %>% 
  mutate(
    unsuitable_area = ca_area - suitable_area_km2, 
    unsuit_pct_ca = 100 - pct_california,
  ) %>% 
  select(
    scenario,
    period,
    suitable_area_km2,
    pct_california,
    unsuitable_area,
    unsuit_pct_ca,
    change_from_baseline_km2,
    pct_change_from_baseline,
    layer
  ) %>% 
  arrange(scenario, period)

table_osmia_area

# Floral resources model results ------------------------------------------

# General stats for floral resources models

# Baseline and future richness patterns



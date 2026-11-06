# libraries
library(sf)
library(terra)
library(furrr)
library(ggridges)
library(tidyterra)
library(tidyverse)
library(patchwork)
library(rnaturalearth)
library(scico)

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

sf_ol <- occs_osmia %>% 
  filter(class != 'background') %>% 
  st_as_sf(coords=c('lon', 'lat'))
st_crs(sf_ol) <- 4326

library(ggspatial)
library(rnaturalearth)

# world <- rnaturalearth::ne_countries(scale = 50)
ggplot() + 
  geom_spatraster(data = base_osmia) + 
  geom_sf(data=sf_ol, colour='gray60',fill='gray10', size=0.8,shape=21, alpha=0.8 )+
  scale_fill_viridis_c(option='H', direction = -1, na.value = NA, alpha = 0.8, 
                       name = 'Suitability') + 
  annotation_scale(style='ticks', location='br') + 
  coord_sf(clip = T)+
  # geom_sf(data=world, colour='gray60',fill=NA, linewidth=0.,alpha=0.8, )+
  theme_minimal()+ 
  theme(panel.background = element_rect(fill=NA, colour='gray80', linewidth=0.5 ))

ggsave(filename = 'figures/F4_general_osmia.png', 
       dpi = 600, width = 23,height = 16, units = 'cm', bg = NULL)

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

ca_ssp585_osmia <- (crop(ssp585_osmia, ca, mask=T)/1000) %>% 
  disagg(fact=8, method='bilinear') %>%
  mask(ca) %>% resample(rout)

raw_change <- round(ca_ssp585_osmia - rout, 3)

pct_change <- raw_change / rout *100



writeRaster(pct_change %>% mask(ca), 
            'outputs/pct_change_high_emission_osmia.tif', overwrite=TRUE)

writeRaster(raw_change %>% mask(ca), 
            'outputs/raw_change_high_emission_osmia.tif', overwrite=TRUE)

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

# full models
fut_rast_full <- rast(list.files(
  "outputs/models/projections/",
  pattern = "\\.tif$",
  full.names = TRUE) %>% str_subset('TSSbinay', negate = T)
  ) 

names(fut_rast_full) <- sources(fut_rast_full) %>% 
  basename() %>% 
  str_remove("_EMmean\\.tif$") %>% 
  str_remove("Osmia.lignaria_") 

fut_rast_full <- fut_rast_full %>% crop(ca, mask=T)/1000

future_table_full <- map_dfr(seq_len(nlyr(fut_rast_full)), \(i) {
  
  layer_name <- names(fut_rast_full)[i]
  
  tibble(
    layer = layer_name,
    suitable_area_km2 = calc_suitable_area(fut_rast_full[[i]], cutoff)
  )
})

# full model results for area
future_table_full <- future_table_full %>% 
  mutate(
    scenario = str_extract(layer, "ssp\\d+|SSP\\d+") %>% str_to_upper(),
    period = case_when(
      str_detect(layer, "2015|2044") ~ "2015--2044",
      str_detect(layer, "2045|2074") ~ "2045--2074",
      str_detect(layer, "2075|2100") ~ "2075--2100",
      TRUE ~ NA_character_
    ), 
    model =str_remove(layer,'_.*')
  ) %>% 
  mutate(
    period = factor(
      period,
      levels = c("2015--2044", "2045--2074", "2075--2100")
    )
  )

# min / max across the individual models, per scenario x period
range_by_model <- future_table_full %>%
  mutate(period = as.character(period)) %>%      # reference table has period as character
  group_by(scenario, period) %>%
  summarise(
    n_models          = n_distinct(model),
    suitable_area_min = min(suitable_area_km2, na.rm = TRUE),
    suitable_area_max = max(suitable_area_km2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # same formula as the reference table, applied to the model extremes
    pct_change_min = 100 * (suitable_area_min - baseline_area) / baseline_area,
    pct_change_max = 100 * (suitable_area_max - baseline_area) / baseline_area,
    across(
      c(suitable_area_min, suitable_area_max, pct_change_min, pct_change_max),
      \(x) round(x, 1)
    )
  )

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

# reference table unchanged; the new columns are only appended
table_osmia_area_range <- table_osmia_area %>%
  left_join(
    range_by_model %>% select(-n_models),
    by = c("scenario", "period")
  ) %>%
  relocate(
    suitable_area_min, suitable_area_max, pct_change_min, pct_change_max,
    .after = pct_change_from_baseline          # keeps `layer` as the last column
  )


# table_osmia_area_range
# table_osmia_area

scen_lab <- c(SSP245 = "Low-emissions\n(SSP245)",
              SSP370 = "Mid-emissions\n(SSP370)",
              SSP585 = "High-emissions\n(SSP585)")
per_lev  <- c("2015-2044", "2045-2074", "2075-2100")
gcms     <- c("CNRM-ESM2-1", "EC-Earth3-Veg", "INM-CM5-0", "MIROC6")

# one row per model x scenario x period, with fixed x positions:
# scenario bin + period offset + small model offset (no random jitter)
pts <- future_table_full %>%
  mutate(period = str_replace(period, '--', '-')) %>% 
  mutate(
    model      = factor(str_extract(layer, str_c(gcms, collapse = "|")), levels = gcms),
    scenario   = factor(scenario, levels = names(scen_lab)),
    period     = factor(as.character(period), levels = per_lev),
    pct_change = 100 * (suitable_area_km2 - baseline_area) / baseline_area,
    x = as.integer(scenario) + (as.integer(period) - 2) * 0.27 +
      (as.integer(model) - 2.5) * 0.045
  )

# ensemble-mean values from the reference table
ens <- table_osmia_area %>%
  mutate(period = str_replace(period, '--', '-')) %>% 
  filter(scenario != "Baseline") %>%
  mutate(
    scenario = factor(scenario, levels = names(scen_lab)),
    period   = factor(as.character(period), levels = per_lev),
    x = as.integer(scenario) + (as.integer(period) - 2) * 0.27
  )

write_csv(pts, 'outputs/future_area_changes_olig.csv')
write_csv(ens, 'outputs/ensemble_future_areas_means_olig.csv')

p <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(data = pts,
             aes(x, pct_change, fill = period, shape = model),
             size = 3.0, colour = NULL, stroke = 0.1) +
  geom_segment(data = ens,
               aes(x = x - 0.1, xend = x + 0.1, colour = period,
                   y = pct_change_from_baseline, yend = pct_change_from_baseline),
               linewidth = 1) +
  scale_x_continuous(breaks = 1:3, labels = scen_lab,
                     expand = expansion(add = 0.5)) +
  scale_fill_scico_d(palette = "batlow", begin = 0.1, alpha = 0.8, end = 0.75, name = "Period") +
  scale_colour_scico_d(palette = "batlow", begin = 0.1, alpha = 0.8, end = 0.75, name = "Period", guide=NULL) +
  
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Global climate model") +
  guides(fill  = guide_legend(order = 1, override.aes = list(shape = 21, size = 3)),
         shape = guide_legend(order = 2, override.aes = list(fill = "grey70", size = 3))) +
  labs(x = NULL, y = "Change in suitable area from baseline (%)") +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.15, 0.25) )

p

# fig. 2 osmia climate suitable area changes
ggsave("figures/fig2_osmia_area_change_models.pdf", p,scale = 1.4,
       width = 12, height = 11, units = "cm", dpi = 300)

# Fig 1C. frequency plots

baseline <- crop(base_osmia, ca, mask=T)
ssp245 <- fut_rast$`ssp245_2075-2100` %>% crop(ca, mask=T)
ssp585 <- fut_rast$`ssp585_2075-2100`%>% crop(ca, mask=T)


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


dats <- 
  tibble(
    baseline = values(baseline, na.rm=T)[,1] %>% unname(), 
 
    ssp245_2100 = values(ssp245, na.rm=T)[,1] %>% unname(),

    ssp585_2100 = values(ssp585, na.rm=T)[,1] %>% unname()
  ) %>% pivot_longer(
    cols = everything()
  ) %>% 
  mutate(Scenario = factor(name, 
                           levels= c(
                             'ssp585_2100',
                             'ssp245_2100',
                             'baseline'
                           )))

ggplot(dats, aes(x = value, y = Scenario, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 1.2, linewidth = 0.1, rel_min_height = 0.005) +
  scale_fill_gradientn(colours = scico(palette = 'nuuk',n = 20, direction = -1))+
  scale_x_continuous(
    breaks = seq(0, 1, by=.2)
  )+
  labs(x = 'Suitability') +
  scale_y_discrete(labels = c('High emissions \n(2075-2100)',
                              'Low emissions \n(2075-2100)',
                              'Baseline'), expand = c(-1.1, 0)) +
  theme_light(base_size = 18) +
  theme(
    legend.position="none"
  )

ggsave('outputs/density_suit_draft.png', height = 14.968, width = 12.205, dpi=600, units = 'cm')  
ggsave('outputs/density_suit_draft.svg', height = 14.968, width = 12.205, dpi=600, units = 'cm')  

# Floral resources model results ------------------------------------------

# General stats for floral resources models
# Load plant list
plant_list <- read_csv('outputs/records/all_species_records_and_native.csv')
natives <- plant_list %>% filter(`CA Native` == 'yes')

# load evaluation metrics
evals <- read_csv("outputs/csv/plants_ensemble_eval.csv") %>% 
  mutate(
    species = str_remove(full.name, "_.*") %>% 
      str_replace_all("\\.", "_"),
    .before = 1
  )

# species modeled
evals %>% pull(species) %>% n_distinct() #out of 71

# total species
natives %>% nrow()

read_csv('outputs/csv/plants_perf_by_PA_size_auc_gated.csv') %>% 
  group_by(PA) %>% 
  summarise(AUC = mean(ROC_median)) #PA4 10,000 AUC 0.864

evals %>% 
  filter(filtered.by == 'ROC') %>% 
  filter(metric.eval=='ROC') %>% 
  pull(calibration) %>% summary()
  
evals %>% 
  filter(filtered.by == 'TSS') %>% 
  filter(metric.eval=='TSS') %>% 
  pull(calibration) %>% summary()

#Load occs
fls_occs <- list.files('outputs/records/plants', '.csv', full.names = T)

#append to database
natives$path <- sapply(natives$species, \(x) fls_occs[grep(str_replace(x, ' ', '_'), fls_occs)], 
                       simplify = T) %>% unname()
natives$sp <- natives$species %>% str_to_lower() %>% str_replace(' ', '_')

# read projection table
csv_mods <- read_csv("outputs/models/plants_projection_table.csv") %>% 
  mutate(
    dir_path = str_remove(path, basename(path)), 
    binary_file = path %>% basename() %>% sub("^([^_]+)_(.+\\.tif)$", "\\1_TSSbinay_\\2", .) %>% str_c(dir_path, .), 
    exists = file.exists(binary_file),
    current_file = paste0("outputs/models/", species, "_NA_cur_ensemble.tif")
  ) %>% 
  left_join(evals, by = "species") %>% 
  left_join(natives %>% select(species= sp, occs_path=path), 
            by = 'species') %>% 
  mutate(occs_path = ifelse(is.na(occs_path), "outputs/records/plants/Arctostaphylos_uva-ursi.csv", 
                            occs_path)) %>% 
  filter(species != 'plagiobothrys_lithocaryus') # endemic outlier

# Baseline and future richness patterns
baseline_plants <- csv_mods %>% 
  distinct(species, current_file, occs_path) 

plan(multisession, workers=18)
plants_bin_10th <- future_map(baseline_plants$species, \(x){
  
  sp_db <- baseline_plants %>% filter(species == x) 

  baseline_map <- rast(sp_db$current_file) %>% mean()/1000

  occs_sp <- read_csv(sp_db$occs_path)
  
  p10_cutoff <- terra::extract(
    baseline_map,
    occs_sp %>% select(decimalLongitude, decimalLatitude)
  ) %>% 
    select(-ID) %>% 
    pull(1) %>% 
    quantile(probs = 0.10, na.rm = TRUE) %>% 
    as.numeric()
    
  ca_plant_baseline <- baseline_map %>% crop(ca, mask=T)
  
  bin_map <- ca_plant_baseline >= p10_cutoff
  # gc()
  return(
    list(
      species = x, 
      cutoff = p10_cutoff,
      ca_baseline = wrap(ca_plant_baseline), 
      ca_binary = wrap(bin_map)
    )
    )  
  })

plan(sequential)


plan(multisession, workers=18)
plants_bin_10th <- future_map(baseline_plants$species, \(x){
  
  sp_db <- baseline_plants %>% filter(species == x) 
  
  baseline_map <- rast(sp_db$current_file) %>% mean()/1000
  
  occs_sp <- read_csv(sp_db$occs_path)
  
  p10_cutoff <- terra::extract(
    baseline_map,
    occs_sp %>% select(decimalLongitude, decimalLatitude)
  ) %>% 
    select(-ID) %>% 
    pull(1) %>% 
    quantile(probs = 0.10, na.rm = TRUE) %>% 
    as.numeric()
  
  ca_plant_baseline <- baseline_map %>% crop(ca, mask=T)
  
  bin_map <- ca_plant_baseline >= p10_cutoff
  # gc()
  return(
    list(
      species = x, 
      cutoff = p10_cutoff,
      ca_baseline = wrap(ca_plant_baseline), 
      ca_binary = wrap(bin_map)
    )
  )  
})

plan(sequential)
# Fig1D California baseline floral resource richness for Osmia lignaria
baseline_floral_rich <- plants_bin_10th %>% map(\(x) unwrap(x$ca_binary) %>% project(ca_ssp585_osmia))

baseline_floral_rich <- baseline_floral_rich %>% rast() %>% sum(na.rm = T)


baseline_floral_rich_resamp <- baseline_floral_rich %>% 
  disagg(fact=8, method='near') %>%
  mask(ca) %>% resample(rout)
# simple representation Fig 1D

ggplot() + 
  geom_spatraster(data = baseline_floral_rich_resamp) + 
  scale_fill_viridis_c(option='D', direction = -1, na.value = NA, alpha = 0.8) + 
  theme_minimal()

#To plot in QGIS
writeRaster(baseline_floral_rich, 'outputs/baseline_floral_richness.tif', overwrite=T)
writeRaster(baseline_floral_rich_resamp, 'outputs/plant_richness_baseline_resamp.tif', overwrite=T)


csv_mods %>% 
  left_join(
    plants_bin_10th %>% map_dfr(\(x) tibble(species = x$species, cutoff = x$cutoff)), by='species'
  ) %>% write_csv('outputs/database_floral_models_10p_th.csv')

# Fig1E percentage change to high-emission scenario
# Baseline and future richness patterns
ssp585_plants <- csv_mods %>% 
  filter(ssp=='ssp585', period == '2075-2100') %>% 
  distinct(species, path) %>% 
  left_join(
    plants_bin_10th %>% map_dfr(\(x) tibble(species = x$species, cutoff = x$cutoff)), by='species'
  )


plan(multisession, workers=18)
ssp_585_plants_bin_10th <- future_map(unique(ssp585_plants$species), \(x){
  # x=ssp585_plants$species[1]
  sp_db <- ssp585_plants %>% filter(species == x) 
  
  high_emissions_map <- (rast(sp_db$path) %>% mean()/1000) %>% crop(ca, mask=T)
  
  ca_high_emissions_bin <- high_emissions_map >= sp_db$cutoff[1]
  wrap(ca_high_emissions_bin)
})

plan(sequential)

ssp585_floral_richn <- map(ssp_585_plants_bin_10th, unwrap) %>% rast() %>% sum(na.rm = T) 

ssp585_floral_richn_resamp <- ssp585_floral_richn %>% 
  disagg(fact=8, method='near') %>%
  mask(ca) %>% resample(rout, method = 'near')

# simple representation Fig 1D
ggplot() + 
  geom_spatraster(data = ssp585_floral_richn_resamp - baseline_floral_rich_resamp) + 
  scale_fill_viridis_c(option='D', direction = -1, na.value = NA, alpha = 0.8) + 
  theme_minimal()

#To plot in QGIS
writeRaster(ssp585_floral_richn, 'outputs/2100_ssp585_floral_floral_richness.tif', overwrite=T)
writeRaster(ssp585_floral_richn_resamp, 'outputs/2100_ssp585_floral_floral_richness_resamp.tif', overwrite=T)

writeRaster(ssp585_floral_richn_resamp - baseline_floral_rich_resamp, 'outputs/plan_ssp585_richness_change_resamp.tif', overwrite=T)



floral_pct_change <- ((ssp585_floral_richn - (baseline_floral_rich+1))/(baseline_floral_rich+1)) * 100

# number of species changes
ggplot() + 
  geom_spatraster(data = (ssp585_floral_richn - baseline_floral_rich)) + 
  scale_fill_gradient2(low='darkred', mid='white', high = 'forestgreen')+
  # scale_fill_viridis_c(option='H', direction = -1, na.value = NA, alpha = 0.8) + 
  theme_minimal()

#percentage changes
ggplot() + 
  geom_spatraster(data = floral_pct_change) + 
  scale_fill_viridis_c(option='H', direction = -1, na.value = NA, alpha = 0.8) + 
  theme_minimal()

writeRaster(floral_pct_change, 'outputs/pct_change_high_emissions_floral_richness.tif', overwrite=T)

# Fig 1F frequency plants

# Baseline and future richness patterns
ssp245_plants <- csv_mods %>% 
  filter(ssp=='ssp245', period == '2075-2100') %>% 
  distinct(species, path) %>% 
  left_join(
    plants_bin_10th %>% map_dfr(\(x) tibble(species = x$species, cutoff = x$cutoff)), by='species'
  )


plan(multisession, workers=18)
ssp_245_plants_bin_10th <- future_map(unique(ssp245_plants$species), \(x){
  # x=ssp245_plants$species[1]
  sp_db <- ssp245_plants %>% filter(species == x) 
  
  high_emissions_map <- (rast(sp_db$path) %>% mean()/1000) %>% crop(ca, mask=T)
  
  ca_high_emissions_bin <- high_emissions_map >= sp_db$cutoff[1]
  wrap(ca_high_emissions_bin)
})

plan(sequential)

ssps24_floral_richn <- map(ssp_245_plants_bin_10th, unwrap) %>% rast() %>% sum(na.rm = T) 

plot(baseline_floral_rich)
plot(ssps24_floral_richn)
plot(ssp585_floral_richn)
dats <- 
  tibble(
    baseline = values(baseline_floral_rich, na.rm=T)[,1] %>% unname(), 

    ssp245_2100 = values(ssps24_floral_richn, na.rm=T)[,1] %>% unname(),

    ssp585_2100 = values(ssp585_floral_richn, na.rm=T)[,1] %>% unname()
  ) %>% pivot_longer(
    cols = everything()
  ) %>% 
  mutate(Scenario = factor(name, 
                           levels= c(
                             'ssp585_2100',
                             'ssp245_2100',
                             'baseline'
                           )))

ggplot(dats, aes(x = value, y = Scenario, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 1, linewidth = 0.1,
                               rel_min_height = 0.005) +
  scale_fill_gradientn(colours = scico(palette = 'imola',n = 20, direction = -1))+
  labs(x = 'Richness') +
  scale_y_discrete(labels = c('High emissions \n(2075-2100)',
                              'Low emissions \n(2075-2100)',
                              'Baseline'), expand = c(-1.1, 0)) +
  theme_light(base_size = 18) +
  theme(
    legend.position="none"
  )

ggsave('outputs/plant_density_suit_draft.png', height = 14.968, width = 12.205, dpi=600, units = 'cm')  
ggsave('outputs/plant_density_suit_draft.svg', height = 14.968, width = 12.205, dpi=600, units = 'cm')  


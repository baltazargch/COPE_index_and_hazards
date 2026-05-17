# libraries
library(sf)
library(terra)
library(furrr)
library(ggridges)
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

range(pct_change)

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
  scale_fill_gradientn(
    colours = c('#d7191c','#fdae61', '#ffffbf', '#abdda4', '#2b83ba') ) +
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

ggsave('outputs/density_suit_draft.png', dpi=600)  
ggsave('outputs/density_suit_draft.svg', dpi=600)  

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
  filter(species != 'plagiobothrys_lithocaryus') # endemic outlier

# Baseline and future richness patterns
plan(multisession, workers=8)
plants_bin_10th <- future_map(csv_mods$species, \(x){
  
  sp_db <- csv_mods %>% filter(species == x)

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
  
  bin_map <- ca_plant_baseline >= cutoff
  
  return(
    list(
      species = x, 
      ca_baseline = ca_plant_baseline, 
      ca_binary = bin_map
    )
    )  
  })

plan(sequential)


# Fig 1F frequency plants
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


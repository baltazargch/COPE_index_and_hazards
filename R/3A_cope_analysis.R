# libraries
library(sf)
library(terra)
library(tidyterra)
library(tidyverse)
library(rnaturalearth)
library(patchwork)

ca <- rnaturalearth::ne_states(country = 'united states of america') %>% 
  filter(name == 'California')

# ---------------------------------------------------------------------------
# INDEX FUNCTIONS
# ---------------------------------------------------------------------------
# Split out wr_index() so omega_r can be computed and plotted on its own,
# not just embedded inside cope_index() after first revision

wr_index <- function(rf, rc, k, v = 0.5, delta = 1, h = 3){
  g_rc <- rc / (rc + h)
  wr <- 1 + v * g_rc * tanh(k * log((rf + delta) / (rc + delta)))
  return(wr)
}

cope_index <- function(sf, rf, rc, k, v = 0.5, delta = 1, h = 3){
  wr <- wr_index(rf, rc, k, v, delta, h)
  cope <- sf * wr
  return(cope)
}

# ---------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------

sf <- rast('outputs/models/projections/means/ssp585_2075-2100_EMmean.tif') %>% 
  crop(ca, mask = TRUE) / 1000

rc <- rast('outputs/baseline_floral_richness.tif')
rf <- rast('outputs/2100_ssp585_floral_floral_richness.tif')

k_val <- 0.5

wr_osmia   <- wr_index(rf, rc, k = k_val, delta = 1, h = 3)
COPE_osmia <- sf * wr_osmia   # equivalent to cope_index(sf, rf, rc, k=k_val, ...)

# ---------------------------------------------------------------------------
# PANEL A: COPE index
# ---------------------------------------------------------------------------

p_cope <- ggplot() +
  geom_spatraster(data = COPE_osmia) +
  scale_fill_terrain_c(name = 'COPE index', na.value = 'transparent') +
  theme_minimal() +
  labs(title = 'A')

writeRaster(COPE_osmia %>%
              disagg(fact = 8, method = 'bilinear') %>%
              mask(ca), 'outputs/COPE_index_osmia.tif', overwrite = TRUE)
writeRaster(share_S_z_masked %>%
              disagg(fact = 8, method = 'bilinear') %>%
              mask(ca), 'outputs/Wr_index_osmia.tif', overwrite = TRUE)

# ---------------------------------------------------------------------------
# PANEL B (v2): standardized spatial attribution
# ---------------------------------------------------------------------------
# Log-additive decomposition, each term's deviation is
# standardized by its own spatial SD before comparing. This removes the
# scale asymmetry from omega_r's inherently further from zero [0.5, 1.5] range vs.
# S_f's close to zero [0-1] range, so "dominance" reflects how unusual a
# pixel is relative to that term's own typical variability, not raw
# log-magnitude.

eps <- 1e-4
S_f_safe <- clamp(sf, lower = eps, values = TRUE)
wr_safe  <- clamp(wr_osmia, lower = eps, values = TRUE)  # wr theoretically >= 0.5, mainly guards edge cases

logS <- log(S_f_safe)
logW <- log(wr_safe)

mean_logS <- global(logS, fun = 'mean', na.rm = TRUE)[1, 1]
mean_logW <- global(logW, fun = 'mean', na.rm = TRUE)[1, 1]

dev_S <- logS - mean_logS
dev_W <- logW - mean_logW

total_dev <- abs(dev_S) + abs(dev_W)
share_S   <- abs(dev_S) / total_dev

sd_logS <- global(logS, fun = 'sd', na.rm = TRUE)[1, 1]
sd_logW <- global(logW, fun = 'sd', na.rm = TRUE)[1, 1]

z_S <- dev_S / sd_logS
z_W <- dev_W / sd_logW

total_z <- abs(z_S) + abs(z_W)
share_S_z <- abs(z_S) / total_z

# Near-average masking, with the standardized total
z_vals      <- values(total_z, na.rm = TRUE)
z_threshold <- quantile(z_vals, probs = 0.05, na.rm = TRUE)
share_S_z_masked <- ifel(total_z < z_threshold, NA, share_S_z)

p_attribution_v2 <- ggplot() +
  geom_spatraster(data = share_S_z_masked) +
  scale_fill_gradient2(
    low = '#8E0152', mid = 'grey90', high = '#276419',
    midpoint = 0.5, limits = c(0, 1),
    name = 'Suitability\nshare (standardized)',
    na.value = 'transparent'
  ) +
  theme_minimal() +
  labs(title = 'B', subtitle = 'Standardized: Red = resource-driven, Green = suitability-driven')

p_attribution_v2

writeRaster(share_S_z_masked %>%
              disagg(fact = 8, method = 'bilinear') %>%
              mask(ca), 'outputs/COPE_attribution_share_Z_osmia.tif', overwrite = TRUE)

# ---------------------------------------------------------------------------
# PANEL C: GCM x SSP uncertainty (SD of COPE across ensemble members)
# ---------------------------------------------------------------------------
# Uses the actual file structure:
#   - Osmia suitability: outputs/models/projections/Osmia.lignaria_{GCM}_{ssp}_{period}_EMmean.tif
#   - Floral richness: uses csv_mods/plants_bin_10th pipeline from 2A_sdm_calc_results.R,
#     generalized over gcm x ssp instead of ssp585 alone

library(furrr)   # for future_map

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

period_target <- '2075-2100'
gcm_list <- c('CNRM-ESM2-1', 'EC-Earth3-Veg', 'INM-CM5-0', 'MIROC6')
ssp_list <- c('ssp245', 'ssp370', 'ssp585')

# -- 1. Generalized richness function (parameterized over ssp and gcm) -----------------

compute_richness_for_scenario <- function(target_ssp, target_gcm, target_period,
                                          csv_mods, plants_bin_10th, ca,
                                          workers = 18) {
  
  scenario_plants <- csv_mods %>%
    filter(ssp == target_ssp, gcm == target_gcm, period == target_period) %>%
    distinct(species, path) %>%
    left_join(
      plants_bin_10th %>% map_dfr(\(x) tibble(species = x$species, cutoff = x$cutoff)),
      by = 'species'
    )
  
  plan(multisession, workers = workers)
  bin_list <- future_map(unique(scenario_plants$species), \(x) {
    sp_db <- scenario_plants %>% filter(species == x)
    r <- (rast(sp_db$path) %>% mean() / 1000) %>% crop(ca, mask = TRUE)
    bin <- r >= sp_db$cutoff[1]
    wrap(bin)
  })
  plan(sequential)
  
  richness <- map(bin_list, unwrap) %>% rast() %>% sum(na.rm = TRUE)
  return(richness)
}

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

# -- 2. Loop over all 12 GCM x SSP combinations, compute COPE for each -------

cope_ensemble_list <- list()

for (gcm in gcm_list) {
  for (ssp in ssp_list) {
    
    combo_name <- paste(gcm, ssp, sep = '_')
    message('Processing: ', combo_name)
    
    # suitability for this GCM x SSP x period
    sf_path <- sprintf('outputs/models/projections/Osmia.lignaria_%s_%s_%s_EMmean.tif',
                       gcm, ssp, period_target)
    sf_gcm <- rast(sf_path) %>% crop(ca, mask = TRUE) / 1000
    
    # floral richness for this GCM x SSP x period
    rf_gcm <- compute_richness_for_scenario(
      target_ssp = ssp, target_gcm = gcm, target_period = period_target,
      csv_mods = csv_mods, plants_bin_10th = plants_bin_10th, ca = ca
    )
    
    # rc (baseline richness) is unchanged across GCM/SSP - loaded once outside this loop
    cope_ensemble_list[[combo_name]] <- cope_index(sf_gcm, rf_gcm, rc, k = k_val, delta = 1, h = 3)
  }
}

cope_ensemble_stack <- rast(cope_ensemble_list)
uncertainty_sd <- app(cope_ensemble_stack, fun = sd, na.rm = TRUE)

p_uncertainty <- ggplot() +
  geom_spatraster(data = uncertainty_sd) +
  scale_fill_viridis_c(option = 'magma', direction = -1, name = 'SD (COPE)',
                       na.value = 'transparent') +
  theme_minimal() +
  labs(title = 'C', subtitle = paste('GCM x SSP spread,', period_target))

writeRaster(uncertainty_sd %>%
              disagg(fact = 8, method = 'bilinear') %>%
              mask(ca), 'outputs/COPE_uncertainty_gcm_ssp.tif', overwrite = TRUE)

# Save the full ensemble stack, useful if you later want the
# variance-decomposition (SSP-dominated vs GCM-dominated) panel
writeRaster(cope_ensemble_stack, 'outputs/COPE_ensemble_gcm_ssp_2075-2100.tif', overwrite = TRUE)

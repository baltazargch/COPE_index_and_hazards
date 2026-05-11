# ============================
# Osmia lignaria – SDM Pipeline
# Tuning → Select best config (AUC gate, TSS rank) → Final models (full data)
# Ensemble from full models → Metrics + Figures
# ============================

# ---- Libraries ----
library(sf)
library(terra)
library(biomod2)
library(tidyverse)
library(fuzzySim)
library(rnaturalearth)
library(rnaturalearthdata)
library(glue)

sf_use_s2(FALSE)

# ---- User paths / I/O ----
dir.create("figures", showWarnings = FALSE)
dir.create("outputs/models", showWarnings = FALSE)
dir.create("outputs/csv", showWarnings = FALSE)

# ---- Parameters ----
set.seed(41546)
auc_min      <- 0.75                # AUC gate for "good discriminant" models
pa_sizes     <- c(2000, 4000, 7000, 9467)  # independent PA sizes
n_cv_reps    <- 4                   # CV replicates
cv_perc      <- 0.80                # train proportion per CV
ncores_small <- 5
ncores_big   <- 14

algos <- c("RF", "GBM", "MAXNET", "GAM", "GLM")  # or MAXNET → MAXENT.Phillips

# ---- Load data ----
# Target group (must include presences & background class column)
tg_os <- readr::read_csv("inputs/records/pres_abs_oslig_target_group.csv",
                         show_col_types = FALSE)

# Presences to sf (EPSG:4326 assumed in CSV)
st_occs <- st_as_sf(tg_os, coords = c("lon", "lat"), crs = 4326)%>%
  filter(class != "background")

# Background points (filter target group == "background")
bg_sf <- st_as_sf(tg_os, coords = c("lon", "lat"), crs = 4326) %>%
  filter(class == "background")

# ---- Env data (historical biovars) ----
# NOTE: if your source file is already EPSG:4326, you can drop rotate/project
hist_fls <- list.files('outputs/NA_bioclim_vars/', '.tif$', full.names = T) %>% 
  stringr::str_subset('historical')

# Calculate mean current conditions
hist_rast <- rast(hist_fls) %>% split(names(.)) %>% map(mean) %>% rast()
names(hist_rast) <- paste('bio', 1:19)
hist_rast <- rotate(hist_rast)
hist_rast <- project(hist_rast, "epsg:4326")

names(hist_rast) <- c(
  "bio1"  = "Ann_Mean_Temp",
  "bio2"  = "Mean_Diurnal_Range",
  "bio3"  = "Isothermality",
  "bio4"  = "Temp_Seasonality",
  "bio5"  = "Max_Temp_Warm_Month",
  "bio6"  = "Min_Temp_Cold_Month",
  "bio7"  = "Ann_Temp_Range", 
  "bio8"  = "Mean_Temp_Wettest_Q", 
  "bio9"  = "Mean_Temp_Driest_Q",
  "bio10" = "Mean_Temp_Warm_Q",
  "bio11" = "Mean_Temp_Cold_Q",
  "bio12" = "Ann_Precip",
  "bio13" = "Precip_Wet_Month",
  "bio14" = "Precip_Dry_Month",
  "bio15" = "Precip_Season",
  "bio16" = "Precip_Wet_Q",
  "bio17" = "Precip_Dry_Q",
  "bio18" = "Precip_Warmest_Q",
  "bio19" = "Precip_Coldest_Q"
)

# ---- Study area: convex buffer of background (100 km in EPSG:5070), then to WGS84 ----
# Using NAD83 / Conus Albers (EPSG:5070) for a planar buffer
bg_aea  <- st_transform(bg_sf, 5070)
sa_poly <- st_union(st_buffer(bg_aea, dist = 100000)) %>% st_transform(4326)

# Clip env to study area
study_mask <- terra::rasterize(terra::vect(sa_poly), hist_rast[[1]], field = 1)
spPreds    <- crop(hist_rast, study_mask, touches = TRUE)

# ---- Variable filtering (collinearity) ----
matPreds <- as.data.frame(spPreds, xy = FALSE, cells = FALSE)
varsCor  <- fuzzySim::corSelect(matPreds, var.cols = 1:ncol(matPreds))
sel_vars <- varsCor$selected.vars
stopifnot(length(sel_vars) >= 2)

# plot(hist_rast)

# ---- Cores depending on area size (km^2) ----
aream <- terra::as.polygons(study_mask, aggregate = TRUE) %>% st_as_sf()
n.cores <- if (as.double(st_area(aream) * 1e-06) <= 20000000) ncores_big else ncores_small

# ---- Build response data (deduplicate to cell centers) ----
# Presence dedup
pres_xy <- st_coordinates(st_occs)
pres_id <- terra::cellFromXY(spPreds[[1]], pres_xy)
pres_xy <- pres_xy[!duplicated(pres_id), , drop = FALSE]

n_pres <- nrow(pres_xy)


# Response vectors: pres = 1, BG = NA (candidates for PA)
resp_xy  <- pres_xy
resp_var <- rep(1, n_pres)


# ---- Format data for biomod2 ----
myBiomodData <- BIOMOD_FormatingData(
  resp.name      = "Osmia_lignaria",
  resp.xy        = resp_xy,
  resp.var       = resp_var,
  expl.var       = spPreds[[sel_vars]],
  PA.nb.rep      = 4,      
  PA.strategy    = "random",
  PA.nb.absences = c(1000,3000,6000,10000),
  filter.raster  = FALSE                     # IMPORTANT: we already deduplicated
)

# ---- Model tuning (with CV) ----
modOL <- BIOMOD_Modeling(
  bm.format     = myBiomodData,
  modeling.id   = "osli.userPA.sized",
  models        = algos,
  CV.strategy   = "random",
  CV.nb.rep     = n_cv_reps,
  CV.perc       = cv_perc,
  OPT.strategy  = "bigboss",
  var.import    = 10,
  metric.eval   = c("ROC", "TSS", "KAPPA"),
  do.full.models= TRUE,                      # build final models on full data
  seed.val      = 123,
  nb.cpu        = n.cores,
  do.progress   = TRUE
)

# =====================================================================
#            TUNING EVALUATIONS → SELECT BEST CONFIG PER ALGO
# =====================================================================

ev_df <- get_evaluations(modOL) %>%
  as_tibble()

# Number of models fitted
length(unique(ev_df$full.name))

# Keep only needed cols and one row per (full.name, metric)
ev_wide <- ev_df %>%
  select(full.name, PA, run, algo, metric.eval, validation) %>%
  tidyr::pivot_wider(
    id_cols    = c(full.name, PA, run, algo),
    names_from = metric.eval,
    values_from= validation,
    values_fn  = max,
    values_fill= NA_real_
  )

biomod2::bm_PlotEvalBoxplot(modOL)

# Gate by AUC (ROC), then summarise TSS per (algo, PA_size)
tuning_by_size <- ev_wide %>%
  filter(!is.na(ROC), ROC >= auc_min) %>%
  group_by(algo, PA) %>%
  summarise(
    TSS_median = median(TSS, na.rm = TRUE),
    ROC_median = median(ROC, na.rm = TRUE),
    KAPPA_median = median(KAPPA, na.rm = TRUE),
    n          = n(),
    .groups    = "drop"
  ) %>%
  mutate(
    PA = case_when(
      PA == 'PA1' ~ 1000,
      PA == 'PA2' ~ 3000,
      PA == 'PA3' ~ 6000,
      PA == 'PA4' ~ 10000, 
      .default =  NA
    )
  ) %>% 
  arrange(PA, algo, desc(TSS_median), desc(ROC_median), desc(KAPPA_median))

# Get some data o the ensemble models
tuning_by_size %>% nrow()
tuning_by_size %>% summary()

readr::write_csv(tuning_by_size, "outputs/csv/tuning_summary_by_size.csv")

# Select the best CONFIG PER ALGO = (PA_size) with max TSS (ties by ROC)
best_size_by_pa <- tuning_by_size %>%
  group_by(algo) %>%
  slice_max(order_by = TSS_median, with_ties = TRUE) %>%
  slice_max(order_by = ROC_median, with_ties = FALSE) %>%
  slice_max(order_by = KAPPA_median, with_ties = FALSE) %>%
  ungroup() %>% 
  arrange(desc(TSS_median), desc(ROC_median), desc(KAPPA_median))

# For reporting too: keep *per-run* records that passed the AUC gate
readr::write_csv(ev_wide %>% filter(!is.na(ROC), ROC >= auc_min),
                 "outputs/csv/tuning_cv_records_auc_gated.csv")

# =====================================================================
#            GET FINAL FULL MODELS FOR THE SELECTED CONFIGS
# =====================================================================

# helper: from a tuning full.name, infer a "base key" (strip RUN#) to match full model
strip_run <- function(x) sub("_RUN\\d+_", "_", x)

built_models <- get_built_models(modOL) %>%
  str_subset("PA1")

# if some algos didn't produce a FULL model name, warn:
if (length(built_models) == 0) {
  stop("No full models detected. Check that BIOMOD_Modeling(..., do.full.models = TRUE) built them.")
}

readr::write_lines(built_models, "outputs/csv/final_full_models_selected.txt")

# =====================================================================
#                     ENSEMBLE FROM FULL MODELS
# =====================================================================

em_algo <- c("EMmean")

bm_ens <- BIOMOD_EnsembleModeling(
  bm.mod        = modOL,
  models.chosen = built_models,
  em.by         = "all",
  em.algo       = em_algo,
  metric.select = "all",
  var.import    = 10,
  seed.val      = 123,
  nb.cpu        = n.cores,
  do.progress   = TRUE
)

saveRDS(bm_ens, file = "outputs/models/osmia_bm_ensemble_model.rds")
saveRDS(modOL, file = "outputs/models/osmia_bm_single_model.rds")

# ---- Ensemble eval summaries (full models only) ----
ens_eval <- get_evaluations(bm_ens) %>% as_tibble()

readr::write_csv(ens_eval, "outputs/csv/ensemble_eval.csv")

# Get some data on ensembles
ens_eval$full.name %>%  unique() %>% length()

# For readability: mean calibration per model/metric
ens_mean <- ens_eval %>%
  group_by(metric.eval, full.name) %>%
  summarise(mean_cal = round(mean(calibration, na.rm = TRUE), 3),
            .groups = "drop") %>%
  select(-full.name) %>% 
  arrange(metric.eval, desc(mean_cal))
readr::write_csv(ens_mean, "outputs/csv/ensemble_eval_means.csv")

# =====================================================================
#                      VARIABLE IMPORTANCE (ensemble models)
# =====================================================================
plot_vi <- biomod2::bm_PlotVarImpBoxplot(
  bm.out  = bm_ens,
  group.by= c("expl.var", "algo", "full.name"),   # per variable x algorithm
  do.plot = FALSE,
  main    = "Variable importance"
)
ens_varimp <- biomod2::bm_PlotVarImpBoxplot(
  bm.out   = bm_ens,
  group.by = c("expl.var", "algo", "full.name"),  # or just "expl.var" if you have one EM method
  do.plot  = FALSE,
  main     = "Ensemble variable importance"
)
ens_varimp$tab %>% 
  group_by(expl.var) %>%
  mutate(mean.val = round(mean(var.imp), 2) *100) %>%
  arrange(desc(mean.val)) %>% 
  select(expl.var, mean.val) %>% 
  distinct()

vi_plot <- plot_vi$plot +
  ggplot2::theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        legend.position = "bottom")
vi_plot
ggsave("figures/fig_sdm_osli_varimp.png", vi_plot,
       width = 8, height = 6, dpi = 300)

# =====================================================================
#                      PROJECTIONS (North America & California)
# =====================================================================

# Current (NA) projection with the same historical stack used for tuning
projOL <- BIOMOD_Projection(
  bm.mod               = modOL,
  models.chosen        = built_models,
  new.env              = spPreds[[sel_vars]],
  proj.name            = "current",
  build.clamping.mask  = TRUE,
  nb.cpu               = n.cores
)

ens_projOL <- BIOMOD_EnsembleForecasting(
  bm.em         = bm_ens,
  bm.proj       = projOL,
  models.chosen = unique(ens_eval$full.name),   # or models_full_selected$full_models
  metric.binary = "TSS",
  nb.cpu        = n.cores
)

# Save NA ensemble raster
ou <- unwrap(ens_projOL@proj.out@val)
terra::writeRaster(ou, "outputs/models/NA_osli_current_ensemble.tif", overwrite = TRUE)

# California stack (provide your own 3-km raster)
models <- c("INM-CM5-0", "EC-Earth3-Veg", "MIROC6", "CNRM-ESM2-1")

periods <- c('1950-2014', '2015-2044', '2045-2074', '2075-2100')

ssp <- c('ssp245', 'ssp370', 'ssp585')

all_bio_files_ca <- list.files('outputs/NA_bioclim_vars', '.tif$', 
                               full.names = T)

pat <- "biovars[_-](\\d{4}[-_]\\d{4})[_-](ssp245|ssp370|ssp585)[_-](.+)\\.tif$"

file_index <- tibble(path = all_bio_files_ca) %>%
  mutate(fname = basename(path)) %>%
  tidyr::extract(fname, into = c("period", "ssp", "gcm"),
                 regex = pat, remove = FALSE) %>%
  mutate(period = str_replace_all(period, "_", "-"),
         gcm = str_remove(gcm, "\\.tif$")) %>%
  filter(!is.na(ssp)) %>%
  arrange(ssp, period, gcm)

# ---- Output dir --------------------------------------------------------------
out_dir <- "outputs/models/projections"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source('R/udf_project_ca_ensemble.R')

proj_tbl <- file_index %>%
  select(path, gcm, ssp, period) %>%
  pmap_dfr(~ project_one(projOL, bm_ens, built_models,unique(ens_eval$full.name), 
                         n.cores, ..1, ..2, ..3, ..4, names(hist_rast), sel_vars, 
                         out_dir))

write_csv(proj_tbl, 'outputs/models/projection_table.csv')

dir.create('outputs/models/projections/means')
dir.create('outputs/models/projections/medians')

# ---- Average across GCMs per (ssp, period) ----------------------------------
mean_tbl <- proj_tbl %>%
  mutate(group = str_c(ssp, '_', period)) %>%
  split(.$group) %>% 
  imap(\(.x, n) {
    
    x <- rast(.x$path) %>% mean()                # stack of EMmean rasters (one per GCM)
    out_mean <- file.path(out_dir, glue("means/{n}_EMmean.tif"))
    writeRaster(x, out_mean, overwrite = TRUE)
    
    x <- rast(.x$path) %>% median()      # mean across GCMs
    out_mean <- file.path(out_dir, glue("medians/{n}_EMmean.tif"))
    writeRaster(x, out_mean, overwrite = TRUE)
    
    tibble(ssp = .x$ssp[1], period = .x$period[1], mean_path = out_mean)
  }) %>% bind_rows()

# write an index CSV of all outputs
write_csv(
  proj_tbl %>% left_join(mean_tbl, by = c("ssp","period")),
  file.path(out_dir, "projection_index.csv")
)
# =====================================================================
#                               MAPS
# =====================================================================

# World basemap via rnaturalearth (land polygons)
wrld <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>% 
  st_crop(spPreds[[1]])

# Quick NA map (pick a layer of the ensemble stack, e.g., 3rd)
p_na <- ggplot() +
  tidyterra::geom_spatraster(data = ou[[3]] / 1000, interpolate = F) +
  geom_sf(data = wrld %>% st_transform(4326), fill = NA, color = "grey60", linewidth = 0.2) +
  scale_fill_viridis_c(option = "H", 
                       name = "Suitability", na.value = NA, 
                       alpha = 0.8, direction = -1) +
  geom_sf(data = st_occs, shape = 20, size = 1.5, stroke = 0.5, 
          color = "grey99", alpha = 0.8) +
  geom_sf(data = st_occs, shape = 20, size = 1.2, stroke = 0.3, 
          color = "grey20", alpha = 0.8) +
  coord_sf(expand = FALSE) +
  ggspatial::annotation_scale(location = "br", 
                              width_hint = 0.25,
                              style = "ticks", line_col="gray20", 
                              text_col = "gray20") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(text = element_text(family = 'sans'),
        panel.background = element_rect(fill='white'),
        legend.position = "right",
        panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.5))

ggsave("figures/fig_sdm_osli_ens_current_NA.png", p_na, width = 8, height = 6, dpi = 400)

# =====================================================================
#              EXTRA: per-PA-size metrics (AUC gate → TSS)
# =====================================================================

# Full per-model CV metrics, gated by AUC, summarised by PA_size & metric
pa_perf <- ev_df %>%
  filter(metric.eval %in% c("ROC", "TSS")) %>%
  select(full.name, PA, run, algo, metric.eval, validation) %>%
  tidyr::pivot_wider(
    id_cols = c(full.name, PA, run, algo),
    names_from = metric.eval, values_from = validation
  ) %>%
  filter(!is.na(ROC), ROC >= auc_min) %>%
  group_by(PA, algo) %>%
  summarise(
    ROC_median = median(ROC, na.rm = TRUE),
    TSS_median = median(TSS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(algo, desc(TSS_median), desc(ROC_median))

readr::write_csv(pa_perf, "outputs/csv/perf_by_PA_size_auc_gated.csv")


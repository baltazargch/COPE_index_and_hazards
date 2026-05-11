# ============================
# Plant resources – SDM Pipeline
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


# UDF ---------------------------------------------------------------------
make_tgb_in_M <- function(
    occ_all, rast_ex, focal, 
    M_focal, n_bg = 2000,
    cap = 1000) {
  
  stopifnot(all(c("species") %in% names(occ_all)))
  M_focal <- st_make_valid(M_focal)
  
  # Target group = other species inside M_focal
  tg <- occ_all %>%
    filter(species != focal) %>%
    group_by(species) %>%
    mutate(.rnd = runif(dplyr::n())) %>% 
    slice_min(.rnd, n = cap, with_ties = FALSE) %>% 
    ungroup() %>%
    st_intersection(M_focal)
  
  rbias <- rast_ex
  
  ll_ras <- rasterize(as.matrix(st_coordinates(tg)), fun='sum',
                      rbias, background=0) 
  
  no_na <- which(!is.na(rbias[]))
  
  ll_rasIDs <- which(values(ll_ras) == 1)
  
  occ_T <- sp::coordinates(raster::raster(ll_ras))[ll_rasIDs,]
  # h <- mean(c(MASS::bandwidth.nrd(occ_T[,2]), MASS::bandwidth.nrd(occ_T[,1])))
  dens <- MASS::kde2d(occ_T[,1], occ_T[,2], h=1.5,
                      n = c(ncol(ll_ras), nrow(ll_ras)), 
                      lims=ext(rbias) %>% as.vector())
  
  biasLayer <- rast(raster::raster(dens)) |> mask(M_focal)
  # plot(biasLayer)
  
  if(!all(dim(biasLayer)[1:2] == dim(rast_ex)[1:2])){
    rbias <- resample(biasLayer, rast_ex, method='bilinear')
  } else {
    rbias[no_na] <- biasLayer[no_na]
  }
  
  bgPoints <- terra::spatSample(rbias, method='weights',
                                n_bg, xy=T, na.rm=T)[, 1:2]
  names(bgPoints) <- c('lon', 'lat')
  # plot(rbias)
  # points(bgPoints, pch='+')
  bgPoints
}

# ---- User paths / I/O ----
dir.create("figures", showWarnings = FALSE)
dir.create("outputs/models", showWarnings = FALSE)
dir.create("outputs/csv", showWarnings = FALSE)

# ---- Parameters ----
set.seed(41546)
auc_min      <- 0.75                # AUC gate for "good discriminant" models
pa_sizes     <- c(2000, 4000, 7000, 10000)  # independent PA sizes
n_cv_reps    <- 4                   # CV replicates
cv_perc      <- 0.80                # train proportion per CV
ncores_small <- 8
ncores_big   <- 16

algos <- c("MAXNET", "GBM", "RF", "GAM", "GLM") #, "GLM",or MAXNET → MAXENT.Phillips

# ---- Load plants occs and list ----
plant_list <- read_csv('outputs/records/all_species_records_and_native.csv')
natives <- plant_list %>% filter(`CA Native` == 'yes')

fls_occs <- list.files('outputs/records/plants', '.csv', full.names = T)

natives$path <- sapply(natives$species, \(x) fls_occs[grep(str_replace(x, ' ', '_'), fls_occs)], 
                       simplify = T) %>% unname()


# Target-group preparation ------------------------------------------------
tg_grp <- map(natives$path, ~read_csv(.x, show_col_types = F)) %>% 
  list_rbind()
tg_grp_sf <- st_as_sf(tg_grp, coords=c('decimalLongitude', 'decimalLatitude'))
st_crs(tg_grp_sf) <- 4326


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

#start per species loop
walk(1:nrow(natives), \(r){
  # r=51
 
  gc()
  occs <- read_csv(natives$path[r], show_col_types = F)
  species <- occs$species[1] %>% janitor::make_clean_names()
  if(file.exists("outputs/csv/plants_perf_by_PA_size_auc_gated.csv")){
    spdone <- read_csv("outputs/csv/plants_perf_by_PA_size_auc_gated.csv", show_col_types = F)
    spdone <- spdone %>% filter(species == species)
    if(any(spdone$species == species)) return()
  }
  
  cat(glue('Analyzing species {species}\n\n'))
  cat(Sys.time() %>% as.character)
  cat('\n\n')
  occs_sf <- st_as_sf(occs, coords=c('decimalLongitude', 'decimalLatitude'))
  st_crs(occs_sf) <- 4326
  
  dups <- terra::extract(hist_rast[[1]], occs_sf, cells=T)$cell %>% duplicated
  
  occs_filt <- occs_sf[!dups, ] %>% st_crop(hist_rast[[1]])
  if(nrow(occs_filt) < 5) return()
  bg_aea  <- st_transform(occs_filt, 5070)
  aream <- st_union(st_buffer(bg_aea, dist = 100000)) %>% st_transform(4326) %>% 
    st_as_sf
  st_crs(aream) <- 4326
  
  envs <- crop(hist_rast, aream, mask=T)
  
  matPreds <- as.data.frame(envs, xy = FALSE, cells = FALSE)
  varsCor  <- fuzzySim::corSelect(matPreds, var.cols = 1:ncol(matPreds))
  sel_vars <- varsCor$selected.vars
  
  stopifnot(length(sel_vars) >= 2)
  envs <-  envs[[ sel_vars ]]
  
  if( sum(!is.na(values(envs[[1]]))) < 15000){
    pa_sizes <- min(sum(!is.na(values(envs[[1]]))) * .8, 10000) %>% floor()
    
    pa_sizes <- c(pa_sizes * 0.2, 
                  pa_sizes * 0.4, 
                  pa_sizes * 0.7, 
                  pa_sizes) %>% floor()
  } else {
    pa_sizes     <- c(2000, 4000, 7000, 10000)  # independent PA sizes
    
  }
  # Example for multiple sizes:
  bg_10000 <- make_tgb_in_M(tg_grp_sf, rast_ex = envs[[1]], 
                            focal = occs_filt$species[1],
                            M_focal = aream, 
                            n_bg = max(pa_sizes))
  
  
  
  # Response vectors: pres = 1, BG = NA (candidates for PA)
  resp_xy  <- occs_filt |> st_coordinates()
  colnames(resp_xy) <- c('lon', 'lat')
  resp_var <- rep(1, nrow(resp_xy))
  
  resp.xy  <- rbind(resp_xy, bg_10000)
  resp.var <- c(rep(1, nrow(occs_filt)), rep(NA, nrow(bg_10000)))
  
  PA.user.table <- matrix(FALSE, nrow = length(resp.var), ncol = length(pa_sizes))
  PA.user.table[1:nrow(occs_filt),] <- TRUE
  bg_idx <- which(is.na(resp.var))
  
  
  set.seed(123)
  for(i in seq_along(pa_sizes)){
    PA.user.table[bg_idx, i][ sample(length(bg_idx), pa_sizes[i]) ]  <- TRUE
  }
  
  
  # ---- Format data for biomod2 ----
  myBiomodData <- BIOMOD_FormatingData(
    resp.name      = janitor::make_clean_names(occs_filt$species[1]),
    resp.xy        = resp.xy,
    resp.var       = resp.var,
    expl.var       = envs[[sel_vars]],
    PA.strategy = "user.defined",
    PA.nb.rep      = 4,      
    PA.user.table = PA.user.table,
    filter.raster  = T                     # IMPORTANT: we already deduplicated
  )
  
  # ---- Model tuning (with CV) ----
  modOL <- BIOMOD_Modeling(
    bm.format     = myBiomodData,
    modeling.id   = paste(occs_filt$species[1], "PA.sized"),
    models        = algos,
    CV.strategy   = "random",
    CV.nb.rep     = n_cv_reps,
    CV.perc       = cv_perc,
    OPT.strategy  = "bigboss",
    var.import    = 10,
    metric.eval   = c("AUCroc", "TSS", "KAPPA"),                    
    seed.val      = 123,
    nb.cpu        = 14,
    do.progress   = F
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
  
  # biomod2::bm_PlotEvalBoxplot(modOL)
  
  # Gate by AUC (AUCroc), then summarise TSS per (algo, PA_size)
  tuning_by_size <- ev_wide %>%
    filter(!is.na(AUCroc), AUCroc >= auc_min) %>%
    group_by(algo, PA) %>%
    summarise(
      TSS_median = median(TSS, na.rm = TRUE),
      ROC_median = median(AUCroc, na.rm = TRUE),
      KAPPA_median = median(KAPPA, na.rm = TRUE),
      n          = n(),
      .groups    = "drop"
    ) %>%
    mutate(
      PA = case_when(
        PA == 'PA1' ~ pa_sizes[1],
        PA == 'PA2' ~ pa_sizes[2],
        PA == 'PA3' ~ pa_sizes[3],
        PA == 'PA4' ~ pa_sizes[4], 
        .default =  NA
      )
    ) %>% 
    arrange(PA, algo, desc(TSS_median), desc(ROC_median), desc(KAPPA_median))
  
  # Get some data o the ensemble models
  tuning_by_size %>% nrow()
  tuning_by_size %>% summary()
  
  readr::write_csv(tuning_by_size, glue("outputs/csv/{species}_tuning_summ_size.csv"))
  
  # Select the best CONFIG PER ALGO = (PA_size) with max TSS (ties by AUCroc)
  best_size_by_pa <- tuning_by_size %>%
    group_by(algo) %>%
    slice_max(order_by = TSS_median, with_ties = TRUE) %>%
    slice_max(order_by = ROC_median, with_ties = FALSE) %>%
    slice_max(order_by = KAPPA_median, with_ties = FALSE) %>%
    ungroup() %>% 
    arrange(desc(TSS_median), desc(ROC_median), desc(KAPPA_median))
  
  # For reporting too: keep *per-run* records that passed the AUC gate
  readr::write_csv(ev_wide %>% filter(!is.na(AUCroc), AUCroc >= auc_min),
                   glue("outputs/csv/{species}_tuning_cv_records_auc_gated.csv"))
  
  # =====================================================================
  #            GET FINAL FULL MODELS FOR THE SELECTED CONFIGS
  # =====================================================================
  
  # helper: from a tuning full.name, infer a "base key" (strip RUN#) to match full model
  strip_run <- function(x) sub("_RUN\\d+_", "_", x)
  
  PA_selected <- tuning_by_size %>% 
    group_by(PA) %>% 
    summarise(meanTSS = mean(TSS_median),
              meanROC = mean(ROC_median),
              meanKAPPA = mean(KAPPA_median),
    ) %>% 
    arrange(meanROC, meanTSS, meanKAPPA) %>% 
    tail(1) %>% pull(PA)
  
  PA_sel <- which(pa_sizes == PA_selected)
  
  built_models <- get_built_models(modOL) %>%
    str_subset(glue("PA{PA_sel}"))
  
  # if some algos didn't produce a FULL model name, warn:
  if (length(built_models) == 0) {
    stop("No full models detected. Check that BIOMOD_Modeling(..., do.full.models = TRUE) built them.")
  }
  
  readr::write_lines(built_models, glue("outputs/csv/{species}_final_full_models_selected.txt"))
  
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
    nb.cpu        = 20,
    do.progress   = F
  )
  
  saveRDS(bm_ens, file = glue("outputs/models/{species}_bm_ensemble_model.rds"))
  saveRDS(modOL, file = glue("outputs/models/{species}_bm_single_model.rds"))
  
  # ---- Ensemble eval summaries (full models only) ----
  ens_eval <- get_evaluations(bm_ens) %>% as_tibble()
  
  readr::write_csv(ens_eval, "outputs/csv/plants_ensemble_eval.csv", append = r!=1)
  
  # Get some data on ensembles
  ens_eval$full.name %>%  unique() %>% length()
  
  # For readability: mean calibration per model/metric
  ens_mean <- ens_eval %>%
    group_by(metric.eval, full.name) %>%
    summarise(mean_cal = round(mean(calibration, na.rm = TRUE), 3),
              .groups = "drop") %>%
    select(-full.name) %>% 
    arrange(metric.eval, desc(mean_cal)) %>% 
    mutate(species = species, .before = metric.eval) %>% 
    distinct() %>% pivot_wider(names_from = metric.eval, values_from = mean_cal)
  
  readr::write_csv(ens_mean, "outputs/csv/ensemble_eval_means.csv", append = r!=1)
  
  # =====================================================================
  #                      VARIABLE IMPORTANCE (ensemble models)
  # =====================================================================
  # bm_ens <- readRDS('outputs/models/acer_macrophyllum_bm_ensemble_model.rds')
  
  ens_varimp <- biomod2::bm_PlotVarImpBoxplot(
    bm.out   = bm_ens,
    group.by = c("expl.var", "algo", "full.name"),  # or just "expl.var" if you have one EM method
    do.plot  = FALSE,
    main     = "Ensemble variable importance"
  )
  
  varimp <- ens_varimp$tab %>% view()
    group_by(expl.var) %>%
    mutate(mean.val = round(mean(var.imp), 2) *100) %>%
    arrange(desc(mean.val)) %>% 
    select(expl.var, mean.val) %>% 
    mutate(species = species, .before =   expl.var ) %>%
    distinct() %>% pivot_wider(names_from = expl.var, values_from = mean.val)
  
  
  varimp[, names(hist_rast)[!names(hist_rast)  %in% colnames(varimp)[-1]]] <- NA
  readr::write_csv(varimp, "outputs/csv/plants_var_importance.csv", append = r!=1)
  
  # =====================================================================
  #                      PROJECTIONS (North America & California)
  # =====================================================================
  
  # Current (NA) projection with the same historical stack used for tuning
  projOL <- BIOMOD_Projection(
    bm.mod               = modOL,
    models.chosen        = built_models,
    new.env              = envs[[sel_vars]],
    proj.name            = "current",
    build.clamping.mask  = TRUE,
    nb.cpu               = 18
  )
  
  ens_projOL <- BIOMOD_EnsembleForecasting(
    bm.em         = bm_ens,
    bm.proj       = projOL,
    models.chosen = unique(ens_eval$full.name),   # or models_full_selected$full_models
    metric.binary = "TSS",
    nb.cpu        = 18
  )
  
  # Save NA ensemble raster
  ou <- unwrap(ens_projOL@proj.out@val)
  terra::writeRaster(ou, glue("outputs/models/{species}_NA_cur_ensemble.tif"), overwrite = TRUE)
  
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
  out_dir <- glue("outputs/models/plants_projections/{species}")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  source('R/udf_project_ca_ensemble.R')
  
  proj_tbl <- file_index %>%
    select(path, gcm, ssp, period) %>%
    pmap_dfr(~ project_one(projOL, bm_ens, built_models,unique(ens_eval$full.name), 
                           20, ..1, ..2, ..3, ..4, names(hist_rast), sel_vars, 
                           out_dir))
  
  
  write_csv(proj_tbl %>% mutate(species=species, .before=gcm), 
            'outputs/models/plants_projection_table.csv', append = r!=1)
  
  dir.create('outputs/models/plants_projections/means')
  dir.create('outputs/models/plants_projections/medians')
  
  # ---- Average across GCMs per (ssp, period) ----------------------------------
  mean_tbl <- proj_tbl %>%
    mutate(group = str_c(ssp, '_', period)) %>%
    split(.$group) %>% 
    imap(\(.x, n) {
      
      x <- rast(.x$path) %>% mean()                # stack of EMmean rasters (one per GCM)
      out_mean <- file.path(glue("outputs/models/plants_projections/means/{species}_{n}_EMmean.tif"))
      writeRaster(x, out_mean, overwrite = TRUE)
      
      x <- rast(.x$path) %>% median()      # mean across GCMs
      out_mean <- file.path(glue("outputs/models/plants_projections/medians/{species}_{n}_EMmean.tif"))
      writeRaster(x, out_mean, overwrite = TRUE)
      
      tibble(ssp = .x$ssp[1], period = .x$period[1], mean_path = out_mean)
    }) %>% bind_rows()
  
  # write an index CSV of all outputs
  write_csv(
    proj_tbl %>% left_join(mean_tbl, by = c("ssp","period")),
    file.path(out_dir, "projection_index.csv")
  )
  
  # =====================================================================
  #              EXTRA: per-PA-size metrics (AUC gate → TSS)
  # =====================================================================
  
  # Full per-model CV metrics, gated by AUC, summarised by PA_size & metric
  pa_perf <- ev_df %>%
    filter(metric.eval %in% c("AUCroc", "TSS")) %>%
    select(full.name, PA, run, algo, metric.eval, validation) %>%
    tidyr::pivot_wider(
      id_cols = c(full.name, PA, run, algo),
      names_from = metric.eval, values_from = validation
    ) %>%
    filter(!is.na(AUCroc), AUCroc >= auc_min) %>%
    group_by(PA, algo) %>%
    summarise(
      ROC_median = median(AUCroc, na.rm = TRUE),
      TSS_median = median(TSS, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(algo, desc(TSS_median), desc(ROC_median)) %>% 
    mutate(species = species, .before = PA) %>% distinct()
  
  readr::write_csv(pa_perf, "outputs/csv/plants_perf_by_PA_size_auc_gated.csv", append = r!=1)
  unlink(bm_ens@sp.name, recursive = T)
  
  # if(r == nrow(natives)){
  #   unlink('rerun.txt')
  # }
  # 
  # if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  #   rstudioapi::restartSession()
  # }
})

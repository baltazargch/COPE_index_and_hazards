library(terra)
library(stars)
library(furrr)
library(tidyverse)
library(tidyterra)
library(rnaturalearth)

fires.files <- list.files('inputs/fires/CNRM-ESM2-1/ssp585/', '.tif', recursive = T, full.names = T)


fires <- rast(fires.files)

time(fires, tstep="years") <- 2001 + 0:99

periods <- list(c(1:14), c(15:44), c(45:74), c(75:100))


ncells <- ncell(fires[[1]])
base = rast(fires[[1]])
vbase <- values(base)

rast_periods <- map(periods, \(pps){
  
  # pps = periods[[1]]
  fires.mat <- fires[[pps]] %>% as.data.frame(na.rm=NA, cells=T) %>% 
    as.matrix() %>% unname()
  
  options(future.globals.maxSize = +Inf)
  plan(multisession, workers=18)
  rle.fires <- future_map(1:nrow(fires.mat), \(x) {
    # x =50000
    rl <- rle(fires.mat[x,-1]) 
    l <- sum(rl$lengths >= 3)
    return(list(
      lenghts = ifelse(l == 0, NA, l),
      intensity = mean(rl$values[ rl$lengths >= 2 ], na.rm=T)
    )
    )
  })
  plan(sequential)
  
  r.length = rle.fires %>% map(~.x$lenghts) %>% list_c() 
  r.intens = rle.fires %>% map(~.x$intensity) %>% list_c()
  
  vbase[!fires.mat[,1]] <- NaN
  vbase[fires.mat[,1]] <- r.length
  duration <- rast(base, vals=vbase)
  vbase[fires.mat[,1]] <- r.intens
  intensity <- rast(base, vals=vbase)
  
  rasts = c(
    intensity, 
    duration
  )
  
  names(rasts) <- c('intenisty', 'duration')
  return(rasts)
})

ca <- ne_states('United States of America') %>% 
  filter(name == 'California') %>% 
  st_transform(crs(rast_periods[[1]][[1]]))

i =4
plot(rast_periods[[i]]$intenisty, col=viridis::magma(100))

plot(ca$geometry, add=T)

dir.create('inputs/fires')
names(rast_periods) <- c('current', 'early', 'mid', 'late')

imap(rast_periods, \(.x, n) writeRaster(.x, paste0('inputs/fires/', n, 'fires_rasters_intens_duration.tif')))



# HEATWAVES  --------------------------------------------------------------
library(seas)
library(furrr)
library(lubridate)
sf_use_s2(FALSE)

## -------------------------
## Helper to flag 4+ day runs in a logical VECTOR
## Returns a logical tibble with events count and maximum duration
## -------------------------
flag_runs_geK <- function(v, k = MIN_RUN) {
  # v = fall_mat[2,]
  # k = MIN_RUN
  
  # if (all(!v) || sum(v, na.rm = TRUE) < k) {
  #   # return(tibble(numberEvents = 0, maxDuration = 0))
  # }
  
  rr <- rle(v)
  
  events <- sum(rr$values & (rr$lengths >= k))
  if(events == 0){
    max <- min <- avg<- 0
  } else {
    max <- max(rr$lengths[ rr$values & (rr$lengths >= k)  ])    
    min <- min(rr$lengths[ rr$values & (rr$lengths >= k)  ])    
    avg <- mean(rr$lengths[ rr$values & (rr$lengths >= k)  ])    
  }
  
  return(c(events, max, min, avg))
  # tibble(numberEvents = events, 
  #             maxDuration = max, 
  #             minDuration = min, 
  #             avgDuration = avg))
}

## ---------------------------------------------------------
## Parameters (match your historical definition)
## ---------------------------------------------------------
THRESH  <- 37     # °C (set 37 if you prefer)
MIN_RUN <- 4      # consecutive days

## ---------------------------------------------------------
## Preload California geometry (you already have this above)
## ---------------------------------------------------------
ne10 <- st_read('/mnt/4TB/GIS/Vectors/NaturalEarth/10m_cultural/ne_10m_admin_1_states_provinces.shp', quiet = TRUE)
cal  <- ne10 %>% filter(name == 'California')

## ---------------------------------------------------------
## Core worker for one scenario directory
## ---------------------------------------------------------
process_loca_dir <- function(dir_path, scenario_label = c("RCP45","RCP85")) {
  scenario_label <- match.arg(scenario_label)
  
  # scenario_label <- 'RCP45'
  # dir_path <- '/mnt/4TB/GIS/Rasters/CalAdapt/LOCA daily max temp/RCP45/'
  
  THRESH  <- 37     # °C (set 37 if you prefer)
  MIN_RUN <- 4      # consecutive days
  
  tif_files <- list.files(dir_path, "\\.tif$", full.names = TRUE)
  stopifnot(length(tif_files) > 0)
  
  # Extract model names from filenames (between 'tasmax_day_' and next '_')
  # Works for paths like .../tasmax_day_MIROC5_r1i1p1_2006-2100_CA.tif, etc.
  get_model <- function(p) gsub(".*tasmax_day_([^_]+).*", "\\1", basename(p))
  models <- vapply(tif_files, get_model, character(1))
  names(tif_files) <- models
  
  # Process each model file
  model_dfs <- vector("list", length(tif_files))
  names(model_dfs) <- models
  
  for (m in seq_along(tif_files)) {
    # m=3
    cat("\nProcessing", scenario_label, "model", m, "of", length(tif_files), ":", names(tif_files)[m], "...\n")
    fl <- tif_files[m]
    
    # We only need metadata first to build date index
    r0    <- rast(fl) - 273.15
    ndays <- nlyr(r0)
    stopifnot(ndays > 0)
    
    # LOCA daily starts at 2006-01-01
    dates <- as.Date("2006-01-01") + 0:(ndays-1)
    
    # Build daily frame and per-year indices
    db_days <- tibble(
      lyridx = 1:ndays,
      date  = dates,
      year  = lubridate::year(dates),
      month = lubridate::month(dates, label = TRUE),
      day   = lubridate::mday(dates)
    )
    
    # mat_hot <- as.matrix(r0 >= THRESH)  # TRUE where Tmax >= THRESH
    
    # Define seasons as indices once
    season_idx <- list(db_days$lyridx
      # spring = db_days$lyridx[db_days$month %in% c('Mar', 'Apr', 'May')]
      # summer = db_days$lyridx[db_days$month %in% c('Jun', 'Jul', 'Aug')],
      # fall   = db_days$lyridx[db_days$month %in% c('Sep', 'Oct', 'Nov')]
    )
    
    hw_fun <- function(x) {
      flag_runs_geK(x >= 37, 4)
    }
    
    terraOptions(memfrac = 0.8)
    
    # plan(multisession, workers = 3)
    # Then per-year × season
    out_dbs <- imap(season_idx, \(x,n){
      # x <- season_idx[[1]]
      xdb <- db_days %>%
        filter(lyridx %in% x) %>%
        group_by(year) %>% dplyr::select(lyridx) %>% 
        split(., .$year)
      
      out_yearly_hw <- imap(xdb, \(z, y){
        # z=xdb[[1]]
        
        rdb <- app(r0[[ z$lyridx ]], hw_fun)
        
        set.names(rdb, c('count', 'max', 'min', 'avg'))
        rdb
        #TODO calculate area per year and area per count
        
        # as.data.frame(rdb) %>% 
        #   as_tibble() %>%
        #   mutate(year = y, .before=1) %>% 
        #   filter(count>0)
      })
      
      out_hw <- map(out_yearly_hw, ~.x$count) %>% rast()
      
      periods <- list(as.character(2006:2014), 
                      as.character(2015:2044), 
                      as.character(2045:2074),
                      as.character(2075:2100))
      
      names(periods) <- c('current', 'early', 'mid', 'late')
      
      hw_periods <- map(periods, \(x){
        # x = periods$late
        
        nms.in <-intersect(as.numeric(names(out_hw)), x)
        
        out <- sum(out_hw[[nms.in]])
        out
      })
      
      hw_periods <- rast(hw_periods)
      hw_periods
      # plot(hw_periods %>% mask(st_transform(ca, 4326)))
    })
    # plan(sequential)
    saveRDS(out_dbs[[1]], paste0('inputs/heatwaves/4DAY_future', 
                            scenario_label, '_',
                            names(tif_files)[m], '_hw.rds'))
  } #END FOR LOOP
} #END FUNCTION

## ---------------------------------------------------------
## Run for RCP 4.5 and 8.5 and write outputs
## ---------------------------------------------------------

plan(multisession, workers = 8)

c('RCP45', 'RCP85') %>% 
  future_walk(., \(rcp){
    library(terra)
    baseDir  <- '/mnt/4TB/GIS/Rasters/CalAdapt/LOCA daily max temp/'
    process_loca_dir(
      dir_path = paste0(baseDir, rcp, '/'),
      scenario_label = rcp
    )
  }, .options = furrr_options(seed=T))

plan(sequential)

ssp245 <- list.files('inputs/heatwaves/','.rds',full.names = T) %>% 
  str_subset('RCP85') %>% map(readRDS)

periods <- c('current', 'early', 'mid', 'late') 

ssp245.periods <- map(periods, ~ rast(map(ssp245, \(y) y[[.x]] )) %>% mean()) %>% 
  rast() %>% crop(st_transform(ca, crs(.)), mask=T)

names(ssp245.periods) <- periods

plot(ssp245.periods)

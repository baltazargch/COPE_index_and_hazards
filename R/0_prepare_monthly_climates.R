library(terra)
library(sf)
library(tidyverse)
library(ncdf4)
library(foreach)
library(doParallel)


# Function to aggregate each model raster to 12 monthly means and convert units
collapse_to_monthly_climatology <- function(r_model, var) {
  # r_model = r_mm
  # Extract month names from layer names
  month_labels <- str_extract(names(r_model), "(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)")
  
  monthly_means <- map(month.abb, function(mon) {
    idx <- which(month_labels == mon)
    # Define seconds in a day to convert units
    seconds_vec <- 86400
    
    if (length(idx) > 0) {
      if(var == 'pr') {
        mean(round(r_model[[idx]] * seconds_vec, 2), na.rm = TRUE)
      } else {
        mean(round(r_model[[idx]],  2), na.rm = TRUE)- 273.15
      }
    } else {
      NA
    }
  })
  
  r_monthly <- rast(monthly_means)
  names(r_monthly) <- month.abb
  return(r_monthly)
}

# Target models and parameters
models= c("INM-CM5-0", "EC-Earth3-Veg", "MIROC6", "CNRM-ESM2-1")

# Time periods and scenarios
periods <- c('1950-2014', 
             '2015-2044', 
             '2045-2074', 
             '2075-2100')
ssp <- c('historical', 'ssp245', 'ssp370', 'ssp585')

# Required variables
vars <- c('pr', 'tasmin', 'tasmax')

# Build case table for each variable × period × scenario
table_cases <- expand_grid(models, vars, periods, ssp) %>% 
  filter((ssp == "historical" & periods == "1950-2014") |
           (ssp != "historical" & periods != "1950-2014")) %>% 
  mutate(case = str_c(vars, '_', periods, '_', ssp))

case_list <- table_cases %>% split(table_cases$case)
# get files names and prepare output directories
clim.files <- list.files('inputs/NAclim/', '.nc$', full.names=T) 


# # get baseline time and models 
# hist <- clim.files[ grep('historical', clim.files) ]

dir.create('outputs/NAmonthly_climates', recursive = TRUE, showWarnings = FALSE)

foreach(i = seq_along(case_list), 
        .packages = c('terra', 'sf', 'tidyverse', 'ncdf4')) %do% {
          # for(i in seq_along(clim.files)) {
          # for(i in 1:2) {
          # i = 2
          # load rasters
          # i=11
          case <- case_list[[i]]
          
          filters <- unique(c(case$vars, case$periods, case$ssp))
          rastfls <- Reduce(function(x, pattern) x[str_detect(x, pattern)], filters, init = clim.files)
          r <- rast(rastfls)
          
          # Determine variable name from filename
          out.name.var <- case$vars %>% unique()
          
          # Get actual year range from raster time dimension
          year.range <- time(r) %>%
            as_date(origin = "1950-01-01") %>%  # If time is in numeric days since 2000-01-01
            year() %>%
            range()
          
          # Combine var + year range (e.g., "tasmax_2000_2020")
          out.name.var <- paste0(out.name.var, "_", str_c(year.range, collapse = "_"), 
                                 '_', case$ssp[1])
          
          if(file.exists(paste0('outputs/NAmonthly_climates/', out.name.var, '.tif'))) {
            cat(paste0(out.name.var, ' already saved\n'))
            return(NULL)
          } else {
            print(paste0('Processing ', out.name.var))
          }
          
          num.years <- diff(year.range)+ 1
          
          stopifnot(nlyr(r) == length(models) * num.years * 12)
          
          names.vec <- rep(models, each = num.years * 12)
          
          times.vec <- time(r) %>%
            month(label = TRUE, abbr = TRUE) %>%
            as.character()
          
          names(r) <- str_c(names.vec, times.vec, sep = '_')
          
          # Split model and month from layer names
          layer_names <- names(r)
          month_names <- str_extract(layer_names, "(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)")
          model_names <- str_remove(layer_names, "_(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$")
          
          # Check structure
          # table(model_names, month_names)
          
          cat('Calculating monthly means\n')
          # Apply to all models
          r.list <- split(r, model_names)
          r.monthly.by.model <- map(r.list, 
                                    ~collapse_to_monthly_climatology(.x,  case$vars %>% unique()))
          
          names(r.monthly.by.model) <- unique(model_names)
          r.monthly.by.model <- r.monthly.by.model %>% imap(\(x,n) {
            names(x) <- str_c(n, '_', names(x)) 
            x
          })
          r_all <- rast(r.monthly.by.model)
          names(r_all) <- str_replace(
            names(r_all),
            "_(\\d{1,2})$",
            \(x) paste0("_", month.abb[as.integer(str_extract(x, "\\d+"))])
          )
          
          cat('Writing means and models\n')
          writeRaster(r_all, paste0('outputs/NAmonthly_climates/', out.name.var, '.tif'))
        }


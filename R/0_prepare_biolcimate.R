library(terra)
library(sf)
library(tidyverse)

# Load your custom biovar calculator using terra
source('R/udf_calculate_biovars.R')  # this should define `terra_biovars()`

#===============================#
# PARAMETERS                    #
#===============================#

# General climate model setup
models <- c("INM-CM5-0", "EC-Earth3-Veg", "MIROC6", "CNRM-ESM2-1")

# Time periods and scenarios
periods <- c('1950-2014', '2015-2044', '2045-2074', '2075-2100')

ssp <- c('historical', 'ssp245', 'ssp370', 'ssp585')

# Required variables
vars <- c('pr', 'tasmin', 'tasmax')

# Input/output directories. Adjusto to NA or CA
input_dir <- 'outputs/NAmonthly_climates'
output_dir <- 'outputs/NA_bioclim_vars'

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#===============================
# PREPARE FILE MATRIX
#===============================

# Build case table for each variable × period × scenario
table_cases <- expand_grid(models, vars, periods, ssp) %>% 
  filter((ssp == "historical" & periods == "1950-2014") |
           (ssp != "historical" & periods != "1950-2014")) %>% 
mutate(ssp = ifelse(periods == '1950-2014', 'historical', ssp)) %>% 
  distinct(.) %>% 
  mutate(case = str_c(periods, '_', ssp)) %>% 
  mutate(fln = str_c(vars, case, sep = '_') %>% str_replace('-', '_') %>% 
           str_c(., '.tif'))

# Split into list of cases: each element corresponds to one time slice + SSP
case_list <- split(table_cases, table_cases$case)

#Load files names
fls <- list.files('outputs/NAmonthly_climates', 
                  '.tif', full.names = T)

#===============================#
# CALCULATE BIOVARS             #
#===============================#
# Loop through each case
bio_list <- map(case_list, \(case_df) {
  message("Processing: ", unique(case_df$case))

  # case_df <- case_list[[1]]
  
  # Read rasters in correct order
  fls <- case_df %>% arrange(match(vars, vars)) %>% pull(fln) %>% unique()

  # IMPORTANT. THIS IS ORDER SPECIFIC. IF FILE NAMES CHANGES, THIS CHANGES.
  # TODO detect variable from filename.
  prec <- rast(str_c(input_dir,'/', fls[1]))
  tmin <- rast(str_c(input_dir,'/', fls[2]))
  tmax <- rast(str_c(input_dir,'/', fls[3]))
  
  # Compute biovars for each model slice
  biovars_by_model <- map(models, \(m) {
    # m <- models[[2]]
    idx <- grep(m, names(prec))
    if (length(idx) != 12) {
      warning("Model ", m, " does not have 12 months; skipping.")
      return(NULL)
    }
    terra_biovars(prec[[idx]], tmin[[idx]], tmax[[idx]])
  })
  
  names(biovars_by_model) <- models
  # plot(biovars_by_model[[4]])
  return(biovars_by_model)
})

#===============================
# WRITE OUTPUTS
#===============================
# Flatten and write
walk2(names(bio_list), bio_list, \(casename, model_stack) {
  walk2(names(model_stack), model_stack, \(model_name, ras) {
    if (!is.null(ras)) {
      outname <- file.path(output_dir, paste0("biovars_", casename, "_", model_name, ".tif"))
      writeRaster(ras, outname, overwrite = TRUE)
    }
  })
})

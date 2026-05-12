#---- Load Required Libraries ----
library(fs)
library(rgbif)
library(terra)
library(tidysdm) 
library(tidyverse)
library(data.table)
library(rnaturalearth)
library(CoordinateCleaner)

c('outputs') %>% map(dir.create)

# Presence data -----------------------------------------------------------

# Osmia lignaria (Megachilidae) 
# GBIF.org (01 October 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.gzs43u
# download the data and put the csv file it in the directory
db <- fread('inputs/records/0037393-240906103802322.csv')
ne_pol <- rnaturalearth::ne_countries(scale = 10)
  
colnames(db)

db_ol <- db %>% filter(species == 'Osmia lignaria') %>% filter(between(year, 1950,2015))
db_ol_cleaned <- clean_coordinates(db_ol, country_ref = ne_pol, 
                                   tests = c("capitals", "centroids", "equal", "gbif", 
                                             "institutions", "outliers", "zeros"))

dir.create('inputs/records', showWarnings = F)

write_csv(db_ol, 'inputs/records/osmia_lignaria_records_august2025.csv')
write_csv(db_ol_cleaned, 'inputs/records/cleaned_osmia_lignaria_records_august2025.csv')

ol_points <- st_as_sf(db_ol_cleaned, coords = c('decimalLongitude', 'decimalLatitude'))

bio <- rast('outputs/NA_bioclim_vars/biovars_1950-2014_historical_CNRM-ESM2-1.tif')[[1]]
bio <-  rotate(bio)
bio <- project(bio, 'epsg:4326')

bbox_am <- ext(-170, -30, -60, 85) 
americas <- ne_pol %>% 
  filter(REGION_UN == 'Americas') %>% 
  st_transform(crs(bio)) %>%  # match CRS
  st_union() %>% 
  st_as_sf() %>% st_crop(bbox_am)


# Create an empty raster template with bio's resolution and americas' extent
template <- rast(
  extent = ext(americas), 
  resolution = res(bio), 
  crs = crs(bio)
)

# Rasterize first
ref <- rasterize(vect(americas), template, field = 1)

osmia <- db_ol_cleaned %>% select(species, basisOfRecord,
                                  coordclean =.summary,
                                  lon=decimalLongitude,
                                  lat=decimalLatitude)
osmia$basisOfRecord %>% table
ext_osm <- terra::extract(ref, osmia %>% st_as_sf(., coords = c('lon', 'lat')), 
                          cells=T)

osmia$inland <- !is.na(ext_osm$layer)
osmia$dups <- duplicated(ext_osm$cell)

osmia_sf <- st_as_sf(osmia, coords=c('lon', 'lat'))
st_crs(osmia_sf) <- 4326

osmia_thin <- thin_by_dist(osmia_sf, dist_min = km2m(5))

osmia_df <- cbind(st_drop_geometry(osmia_sf), st_coordinates(osmia_sf))
osmia_df %>% 
  filter(inland) %>% 
  select(
    species:dups, 
    lon=X, lat=Y
  ) %>% write_csv(., 'inputs/records/clean_thin_osmia_lignaria.csv')

# background data target group approach
olig <- read_csv('inputs/records/clean_thin_osmia_lignaria.csv') 
olig <- olig %>% filter(coordclean, inland, !dups)

olig_pp <- st_as_sf(olig, coords= c('lon', 'lat'))
st_crs(olig_pp) <- 4326

osmia <- db %>% filter(genus=='Osmia')
osmia_pp <- osmia %>% st_as_sf(., coords = c('decimalLongitude', 'decimalLatitude'))
st_crs(osmia_pp) <- 4326

# which regions are in the data?
v <- refv <- values(ref)
v[ !is.na(v) ] <- 1
values(ref) <- v
plot(ref)
points(olig[, c('lon', 'lat')])

osmia_bg <- rasterize(osmia_pp, terra::aggregate(ref, 20), fun = "count")
osmia_bg <- terra::disagg(osmia_bg, 20) %>% crop(bio)
plot(osmia_bg)

set.seed(1234567)
olig_w_bg <- sample_background(data = olig_pp, 
                               raster = osmia_bg,
                               n = 10000,
                               method = "bias",
                               class_label = "background",
                               return_pres = TRUE)

st_crs(olig_w_bg) <- 4326

inland  <- terra::extract(bio, olig_w_bg, cells=T)
olig_w_bg$inland <- !is.na(inland$bio1)

olig_w_bg <- olig_w_bg %>% filter(inland)

records_pseudo <- cbind(st_drop_geometry(olig_w_bg), 
                        st_coordinates(olig_w_bg)) %>% 
  dplyr::rename(
    lon=X, lat=Y
  )

write_csv(records_pseudo, 'inputs/records/pres_abs_oslig_target_group.csv')

# Floral resources (California natives) 
# GBIF.org (05 August 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.9c7p3p 
# GBIF.org (05 August 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.5pj7bq

#---- Script: Batch GBIF Downloads with DOI Tracking and Per-Species Outputs ----
# Description:
# Downloads GBIF occurrences using `occ_download()` in batches (default 50 species),
# respecting the limit of 3 concurrent pending downloads. Tracks DOIs for reproducibility,
# cleans coordinates, and saves cleaned per-species CSVs for SDM applications.
#
# Inputs:
# - inputs/plant_species.csv and inputs/animal_species.csv (with `prob_valid`)
#
# Outputs:
# - Cleaned occurrence CSVs per species under: outputs/records/plants/ and .../animals/
# - Metadata log with binomial, usageKey, downloadKey, DOI: outputs/gbif_download_log.csv
#
# Requirements:
# - GBIF credentials stored in .Renviron (GBIF_USER, GBIF_PWD, GBIF_EMAIL)
# - Packages: tidyverse, rgbif, CoordinateCleaner, fs


# List of plants
plants_osmia_resources <- read_delim("inputs/plants_osmia_resources.csv", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
plants_osmia_resources %>% 
  filter(`CA Native` == 'yes') %>% 
  count


#---- Parameters ----
batch_size <- 100
output_log <- "outputs/gbif_download_log.csv"
max_pending <- 3

#---- Create Output Directories ----
dir_create("outputs/records/raw")
dir_create("outputs/records/plants")
dir_create("outputs/tmp_downloads")

#---- Utility Functions ----

##---- Robust status check with retry loop ----
check_status_with_retry <- function(key, max_attempts = 10, wait_secs = 30) {
  attempt <- 1
  while (attempt <= max_attempts) {
    result <- tryCatch({
      meta <- occ_download_meta(key)
      return(meta$status)
    }, error = function(e) {
      message(glue::glue("⚠️ Error on attempt {attempt} for key {key}: {e$message}"))
      return(NA)
    })
    
    if (!is.na(result)) return(result)
    
    Sys.sleep(wait_secs)
    attempt <- attempt + 1
  }
  return("UNKNOWN")
}


##---- Match names to GBIF backbone ----
match_names <- function(species_vec) {
  name_backbone_checklist(name = species_vec) %>%
    filter(!is.na(usageKey)) %>%
    select(original = canonicalName, species, status, usageKey)
}

##---- Submit a batch download ----
submit_download_batch <- function(binomials, usageKeys, user, pwd, email) {
  download_key <- occ_download(
    pred_in("taxonKey", usageKeys),
    pred("hasCoordinate", TRUE),
    format = "SIMPLE_CSV"
  )
  tibble(
    binomial = binomials,
    usageKey = usageKeys,
    downloadKey = download_key,
    date_requested = Sys.time()
  )
}

##---- Clean and split downloaded file by species ----
process_download <- function(downloadKey, group_df) {
  f <- occ_download_get(downloadKey, overwrite = TRUE) %>%
    occ_download_import()
  # group_df <- submission_log[[1]]
  
  f %>% select(species, decimalLatitude, decimalLongitude, countryCode,
               uncertainty=coordinateUncertaintyInMeters,
               year, basisOfRecord, gbifID, taxonKey) %>% 
    write_csv(glue::glue('outputs/records/raw/{downloadKey}_raw.csv'))
  
  species_list <- group_df$species
  
  cleaned <- f %>%
    filter(species %in% species_list,
           !is.na(decimalLatitude),
           !is.na(decimalLongitude)) %>%
    clean_coordinates(
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      species = "species",
      test=c('seas', 'institutions', 'gbif', 'equal','duplicates'),
      value = "clean"
    ) %>%
    select(species, decimalLatitude, decimalLongitude, countryCode,
           uncertainty=coordinateUncertaintyInMeters,
           year, basisOfRecord, gbifID, taxonKey) %>%
    rename(usageKey = taxonKey)
  
  walk(unique(cleaned$species), function(sp) {
    df <- cleaned %>% filter(species == sp)
    folder <- if (sp %in% group_df$species[group_df$group == "animals"]) "animals" else "plants"
    fname <- str_replace_all(sp, " ", "_")
    write_csv(df, file = glue::glue("outputs/records/{folder}/{fname}.csv"))
  })
}

#---- Load and Prepare Species Lists ----
plants <- plants_osmia_resources$Plant

plants_match <- match_names(plants) %>% mutate(group = "plants")

all_species <- plants_match %>% 
  filter(!is.na(original))

# notin <- all_species[-match(sp_done, all_species$original),]

#---- Split into Batches ----
batches <- split(all_species, ceiling(seq_along(all_species$original) / batch_size))

#---- Sequential Submission with Pending Limit ----
submission_log <- list()

for (i in seq_along(batches)) {
  
  # Submit batch
  batch <- batches[[i]]
  cat(glue::glue("Submitting batch {i}/{length(batches)} with {nrow(batch)} species...\n\n"))
  
  dl <- submit_download_batch(batch$species, batch$usageKey)
  
  completed <- FALSE
  
  while (!completed) {
    Sys.sleep(120)
    meta <- occ_download_meta(dl$downloadKey[1])
    statuses <- meta$status
    doi <- meta$doi
    completed <- statuses == "SUCCEEDED"
    cat(glue::glue("[{Sys.time()}] Completed: {sum(completed)} / 1\n\n"))
  }
  
  dl <- dl %>%  rename(species = binomial)
  submission_log[[i]] <- left_join(batch, dl, by = c("species", "usageKey"))
  submission_log[[i]]$doi <- doi
  
  message("Saving current download...")
  
  #---- Process and Clean All Finished Downloads ----
  process_download(dl$downloadKey[1], submission_log[[i]])
  
  species_out <- setdiff(batch$species, 
                         list.files('outputs/records/', '.csv', recursive = T) %>% 
                           basename() %>% 
                           str_remove('.csv') %>% 
                           str_trim() %>% 
                           str_squish() %>% str_replace('_', ' '))
  
  write_lines(species_out, file = 'outputs/records/log_species_not_found.txt', 
              append = i != 1)
  
}

#---- Save Metadata Log ----



log_df <- do.call(rbind, submission_log)

write_csv(log_df, output_log)

message("All GBIF occurrences downloaded, cleaned, and saved per species.")

downloads <- occ_download_list(limit=100)

raws <- list.files('outputs/records/raw/', '.csv') %>% str_remove('_raw.csv')

download_data <- downloads$results %>% filter(key %in% raws)

download_data %>% select(key, doi, status, downloadLink, 
                         created, totalRecords, request.format) %>% 
  write_csv('outputs/records/metadata_downloads_dois.csv')

sp_in <- list.files('outputs/records/', '.csv', recursive = TRUE)[-1] %>% 
  enframe(name = NULL, value = "path") %>%
  separate_wider_delim(path,
                       names = c('group', 'file'), delim = "/", 
                       cols_remove = FALSE) %>%
  filter(group != 'raw') %>% 
  mutate(
    especie = file %>% str_remove("\\.csv$") %>% str_replace("_", " "),
    path = str_c('outputs/records/', path)
  ) %>%
  select(group, especie, path) %>%
  mutate(
    nobs = map_int(path, ~ {
      read_csv(.x, show_col_types = F, progress = F) %>%
        distinct(decimalLatitude, decimalLongitude) %>%
        nrow()
    })
  )

sp_in <- sp_in %>% filter(group != 'raw')

sp_in %>% filter(nobs < 6) %>% count(group)

sp_all <- c(plants) %>% 
  match_names() %>% 
  filter(!is.na(original), !is.na(species)) 

sp_all <- sp_all %>% 
  distinct(species, .keep_all = T) %>% 
  left_join(plants_osmia_resources, by=c('original' = "Plant"))

write_csv(sp_all, 'outputs/records/all_species_records_and_native.csv')

sp_all$species[!sp_all$species %in% sp_in$especie] %>% 
  write_lines('outputs/records/species_not_in_gbif.txt')


library(data.table)
library(tidyverse)
library(CoordinateCleaner)
library(sf)
library(scico) 
library(tidysdm)
library(terra)
library(tidyterra)
library(ggspatial)

sf_use_s2(F)


db <- fread('/mnt/4TB/Databases/GBIF/Megachilidae/0037393-240906103802322.csv')
ne_pol <- read_sf('/mnt/4TB/GIS/Vectors/NaturalEarth/10m_cultural/ne_10m_admin_0_countries.shp')
colnames(db)

db_ol <- db %>% filter(species == 'Osmia lignaria') %>% filter(between(year, 1950,2015))
db_ol_cleaned <- clean_coordinates(db_ol, country_ref = ne_pol, 
                                   tests = c("capitals", "centroids", "equal", "gbif", 
                                             "institutions", "outliers", "zeros"))

dir.create('data/records', showWarnings = F)

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
                               return_pres = TRUE) #%>% 
  # nest_by(class) %>%
  # mutate(data = list(if (class == "background") {
  #   sample_frac(data, 0.7)
  # } else {
  #   data
  # })) %>%
  # unnest(data) %>% st_as_sf(sf_column_name='geometry')

st_crs(olig_w_bg) <- 4326

inland  <- terra::extract(bio, olig_w_bg, cells=T)
olig_w_bg$inland <- !is.na(inland$bio1)

olig_w_bg <- olig_w_bg %>% filter(inland)

ppref_america <- project(ref, "EPSG:4269")
ppolig_w_bg <- st_transform(olig_w_bg, crs = 4269)


ggplot() +
  # Light background raster to show extent (e.g., elevation, mask, etc.)
  geom_spatraster(data = ppref_america, aes(fill = layer), alpha = 1, show.legend = F) +
  scale_fill_wiki_c(na.value = 'transparent')+
  # Background points
  geom_sf(data = ppolig_w_bg %>% filter(class == "background"),
          aes(color = "Background"),
          size = 0.5, alpha = 0.5) +
  
  # Presence points
  geom_sf(data = ppolig_w_bg %>% filter(class != "background"),
          aes(color = "Presence"),
          size = 0.4, alpha = 0.9) +
  
  # Scale bar
  annotation_scale(style='ticks', pad_x = unit(0.3, "cm"),
                   pad_y = unit(0.3, "cm"),)+
  
  # Coordinate system
  coord_sf(crs = "EPSG:4269") +
  
  # Legend and colors
  scale_color_manual(
    name = NULL,
    labels = c('Pseudo-absence', 'Presence'),
    values = c("Presence" = "darkgreen", "Background" = "gray40")
  ) +
  
  # Theme styling
  # Theme and embedded legend
  theme_minimal(base_size = 11) +
  theme(
    panel.border = element_rect(fill = NA, color = "gray30"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    
    # ✅ NEW syntax for inside legend placement
    legend.position = 'inside',
    legend.position.inside = c(0.98, 0.02),
    legend.justification = c(1, 0),
    
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.text = element_text(size = 10)
  )


records_pseudo <- cbind(st_drop_geometry(olig_w_bg), 
                        st_coordinates(olig_w_bg)) %>% 
  dplyr::rename(
    lon=X, lat=Y
  )

write_csv(records_pseudo, 'inputs/records/pres_abs_oslig_target_group.csv')

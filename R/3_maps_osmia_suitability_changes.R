library(terra)
library(sf)
library(tidyverse)

# ---- Load data ----
# Target group (must include presences & background class column)
tg_os <- readr::read_csv("inputs/records/pres_abs_oslig_target_group.csv",
                         show_col_types = FALSE)

# Presences to sf (EPSG:4326 assumed in CSV)
st_occs <- st_as_sf(tg_os, coords = c("lon", "lat"), crs = 4326)%>%
  filter(class != "background")

r <-  rast("outputs/models/NA_osli_current_ensemble.tif")

flsmodels <- list.files('outputs/models/projections/', '.tif$', full.names = T)

periods <- c('2015-2044', 
             '2045-2074', 
             '2075-2100')

ssp <- c('ssp245', 'ssp370', 'ssp585')

dbperiods <- expand_grid(periods=periods, ssp=ssp)

for(i in 1:nrow(dbperiods)){
  row <- dbperiods[i,]
  
  flsin <- flsmodels[grep(row$periods, flsmodels)]
  flsin <- flsin[grep(row$ssp,flsin)]
  
  r1 <- rast(flsin)
  
  # World basemap via rnaturalearth (land polygons)
  wrld_org <- rnaturalearth::ne_states(country = 'united states of america', 
                                       returnclass = "sf")
  
  ca <- wrld_org %>% filter(name_en == 'California')
  
  # World basemap via rnaturalearth (land polygons)
  wrld <- wrld_org %>% 
    st_crop(extend(ext(crop(r1 - r[[1]], ca)), 0.5))
  
  difr <- crop(r1 - r[[1]], ca, mask=T) / 1000
  
  difr <- ifel(difr < -0, -1, ifel(difr > 0, 1, 0))
  # Small helper for consistent breaks/labels across legend(s)
  cat_breaks <- c("-1", "1", "0")
  cat_labels <- c("Loss", "Gain", "Stable")
  
  # 1) Define order and relabel the categorical raster correctly (note the 'ID' col)
  difr <- modal(difr)
  difr <- as.factor(difr)
  levels(difr) <- data.frame(
    ID    = c(-1, 1, 0),
    class = cat_labels
  )
  
  
  library(patchwork)
  
  # --- (A) Area by class (km²) -----------------------------------------------
  # difr: SpatRaster with values -1 (loss), 0 (stable), 1 (gain), already cropped to CA
  # compute per-cell area (km²) respecting latitude
  cell_km2 <- terra::cellSize(difr, unit = "km")           # SpatRaster of cell areas
  cell_km2 <- terra::mask(cell_km2, difr)                  # align NA to difr
  
  area_df <- tibble(
    class_raw = as.vector(terra::values(difr)),
    area_km2  = as.numeric(terra::values(cell_km2))
  ) %>%
    filter(!is.na(class_raw), !is.na(area_km2)) %>%
    mutate(class = factor(
      dplyr::recode(class_raw, `-1` = "Loss", `0` = "Stable", `1` = "Gain"),
      levels = c("Loss", "Gain", "Stable")
    )) %>%
    group_by(class) %>%
    summarise(area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop") %>%
    mutate(pct = 100 * area_km2 / sum(area_km2)) %>% 
    arrange(pct)
  
  
  pal <- setNames(c("#5B1A18",  # Loss  (orange)
                    "#0B775E",  # Gain  (blue)
                    "darkorange2"), cat_labels)
  
  # Set all to desired CSR
  difrp <- disagg(difr, fact=8) %>% mask(ca) %>% project(crs('epsg:4269')) 
  wrldp <- st_transform(wrld, crs=4269)
  
  # --- (B) Main map (your plot, lightly tidied) -------------------------------
  p_map <- ggplot() +
    geom_sf(
      data = wrldp,colour='gray50', fill='gray88',
      show.legend = F
    ) +
    tidyterra::geom_spatraster(data = difrp, interpolate = FALSE,
                               alpha=0.8, na.rm = T) +
    scale_fill_manual(
      values = pal,
      limits = cat_labels, breaks = cat_labels,
      name = "Suitability change",
      na.translate = FALSE   # <- removes the NA legend key
    ) + 
    coord_sf(expand = F, xlim = c(-124.5, -114), ylim = c(32.5, 42.2), crs = 4269) +
    ggspatial::annotation_scale(text_cex=1,
                                location = "bl", width_hint = 0.25, style = "ticks",
                                line_col = "gray20", text_col = "gray20"
    ) +
    labs(x = NULL, y = NULL, title = paste0(row$periods, ' ', row$ssp)) +
    theme_minimal(base_size = 16) +
    theme(
      text = element_text(family = "sans"),
      panel.background = element_rect(fill = "#D9E8F5"),
      legend.position = "bottom",
      panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.5), 
    )
  
  # --- (C) Inset barplot ------------------------------------------------------
  p_bar <- ggplot(area_df, aes(x = class, y = area_km2 / 1000, fill = class)) +
    geom_col(width = 0.75, show.legend = F) +
    geom_text(aes(label = sprintf("%.0f%%", pct)),
              vjust = -0.35, size = 3) +
    scale_fill_manual(
      values = pal,
      name = "Suitability change",
      na.translate = FALSE
    ) + 
    labs(x = NULL, y = expression(Area~(km^2)~x~10^3)) +
    theme_light(base_size = 12) +
    scale_y_continuous(expand = expansion(mult = c(0., 0.15))) +
    coord_cartesian(clip = "off")+
    theme(
      plot.margin = margin(4, 4, 4, 4),
      axis.title.y = element_text(margin = margin(r = 2)),
      # axis.text.x = element_text(angle = 0, vjust = 1),
      axis.text.x = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(),
      plot.background = element_rect(color='gray30', fill = 'white')
    )
  
  # --- (D) Compose with patchwork (top-right inset) ---------------------------
  # Adjust the numbers [0..1] to fine-tune placement if needed.
  final_plot <- p_map + patchwork::inset_element(
    p_bar,
    left   = 0.63,  # x0
    bottom = 0.62,  # y0
    right  = 0.99,  # x1
    top    = 0.98,  # y1
    align_to = "plot"
  ) 
  
  # Save
  ggsave(plot = final_plot, 
         filename = paste0('figures/suit_change_',row$periods, '_', row$ssp, '.png'), 
         dpi=600,
         width = 14, height = 16, scale=1.5, units = 'cm'
  )
  dir.create('figures/svg_suit_change', showWarnings = F)
  # Save
  ggsave(plot = final_plot, 
         filename = paste0('figures/svg_suit_change/', 
                           row$periods, '_', row$ssp, '.svg'), 
         dpi=600,
         width = 14, height = 16, scale=1.2, units = 'cm'
  )
}

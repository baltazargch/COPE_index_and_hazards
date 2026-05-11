library(terra)
library(sf)
library(tidyverse)

# World basemap via rnaturalearth (land polygons)
wrld_org <- rnaturalearth::ne_states(country = 'united states of america', 
                                     returnclass = "sf")

ca <- wrld_org %>% filter(name_en == 'California')

#base layer to extend when plants are not in whole California
r <-  rast("outputs/models/NA_osli_current_ensemble.tif") %>% 
  crop(ca, mask=T) / 1000
wrld <- wrld_org %>% 
  st_crop(extend(ext(r), 0.5))

res <- read_csv("outputs/final_maxent/csv/final_models_metrics.csv")

all_fls <- list.files('outputs/final_maxent/rasters/', 'p10.tif$', 
                      full.names = T)

current_richness <- all_fls[ !str_detect(all_fls, 'ssp') ]


c_rich <- map(current_richness, ~ rast(.x) %>% 
                crop(ca, mask=T) %>% 
                resample(r, method='near')) 

c_rich <- c_rich %>% rast() %>% sum(na.rm = T)
plot(c_rich)
writeRaster(disagg(c_rich, fact=8, method='near') %>% mask(ca), 'figures/plants_current_richness.tif', 
            overwrite=T)

cases <- colnames(res)[13:ncol(res)]

rastOut <- list()
rastRichFut <- list()
for(c in seq_along(cases)){
  # c = 1
  fut_richness <- all_fls[ str_detect(all_fls, cases[c]) ]
  
  
  f_rich <- map(fut_richness, ~ rast(.x) %>% 
                  resample(r, method='near')) 
  
  f_rich <- f_rich %>% rast() %>% sum(na.rm = T) 
  
  rastRichFut[[c]] <- disagg(f_rich, fact=8, method='near') %>% mask(ca)
  
  d_rich <- (f_rich - c_rich) / c_rich * 100 
  cur0 <- c_rich == 0
  
  d_rich[cur0 & (f_rich == 0)] <- 0
  names(d_rich) <- 'Change%'
  
  library(patchwork)
  
  # Set all to desired CSR
  difrp <- disagg(d_rich, # clamp(d_rich, lower = -100, upper = 100), 
                  fact=8, method='bilinear') %>% 
    mask(ca) %>% project(crs('epsg:4269')) 
  rastOut[[c]] <- difrp
  
  # wrldp <- st_transform(wrld, crs=4269)
  # 
  # # make the color range 
  # # Keep your actual range (no symmetry), and label ticks as percentages
  # rng   <- as.numeric(terra::global(difrp, fun = "range", na.rm = TRUE))   # c(min, max)
  # brks  <- pretty(rng, n = 6)
  # 
  # # --- (B) Main map (your plot, lightly tidied) -------------------------------
  # p_map <- ggplot() +
  #   geom_sf(data = wrldp, colour = "gray50", fill = "gray88", show.legend = FALSE) +
  #   tidyterra::geom_spatraster(data = difrp, 
  #                              interpolate = T, alpha = 0.9, na.rm = TRUE) +
  #   scale_fill_gradientn(
  #     colours  = c('darkred', "#B8574E", "#F2F2F2", "#228833", 'darkgreen'),
  #     # limits = c(-100, 100),
  #     labels = scales::label_number(suffix = "%"), name = "Richness change",
  #     oob = scales::squish, na.value = NA
  #   )+
  #   coord_sf(expand = F, xlim = c(-124.5, -114), ylim = c(32.5, 42.2), crs = 4269) +
  #   ggspatial::annotation_scale(text_cex=1,
  #                               location = "bl", width_hint = 0.25, style = "ticks",
  #                               line_col = "gray20", text_col = "gray20"
  #   ) +
  #   labs(x = NULL, y = NULL, 
  #        title = cases[c]) +
  #   theme_minimal(base_size = 16) +
  #   theme(
  #     text = element_text(family = "sans"),
  #     panel.background = element_rect(fill = "#D9E8F5"),
  #     legend.position = "right",
  #     panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.5), 
  #   )
  # 
  # p_map
  # # --- (C) Inset barplot ------------------------------------------------------
  # # area-weighted bins (km²)
  # cell_km2 <- terra::mask(terra::cellSize(difrp, unit="km"), difrp)
  # df <- tibble(change = as.numeric(terra::values(difrp)),
  #              area = as.numeric(terra::values(cell_km2))) |>
  #   filter(is.finite(change), is.finite(area))
  # 
  # breaks <- c(-100, -50, -10, 10, 50, 100)
  # labs   <- c("≤ -50%", "-50 to -10%", "-10 to 10%", "10 to 50%", "≥ 50%")
  # 
  # bin_df <- df |>
  #   mutate(bin = cut(change, breaks = breaks, labels = labs, include.lowest = TRUE)) |>
  #   count(bin, wt = area, name = "area_km2") |>
  #   mutate(pct = 100 * area_km2 / sum(area_km2),
  #          bin = factor(bin, levels = labs))
  # 
  # p_inset <- ggplot(bin_df, aes(x = bin, y = pct, fill = bin)) +
  #   geom_col(width = 0.8, show.legend = FALSE) +
  #   geom_text(aes(label = sprintf("%.0f%%", pct)), vjust = -0.25, size = 3) +
  #   scale_fill_manual(values = c("#B8574E","#E6B0A4","#F2F2F2","#A7C7B7","#228833")) +
  #   scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
  #   labs(x = NULL, y = "Area (%)") +
  #   theme_light(base_size = 10) +
  #   theme(
  #     plot.margin = margin(4, 4, 4, 4),
  #     axis.title.y = element_text(margin = margin(r = 2)),
  #     # axis.text.x = element_text(angle = 0, vjust = 1),
  #     panel.grid.minor = element_blank(), 
  #     panel.grid.major = element_blank(),
  #     plot.background = element_rect(color='gray30', fill = 'white'),
  #     axis.text.x = element_text(angle = 35, hjust = 1)
  #     )
  # 
  # 
  # # --- (D) Compose with patchwork (top-right inset) ---------------------------
  # # Adjust the numbers [0..1] to fine-tune placement if needed.
  # final_plot <- p_map + patchwork::inset_element(
  #   p_inset,
  #   left   = 0.63,  # x0
  #   bottom = 0.62,  # y0
  #   right  = 0.99,  # x1
  #   top    = 0.98,  # y1
  #   align_to = "plot"
  # ) 
  # 
  # # Save
  # ggsave(plot = final_plot, 
  #        filename = paste0('figures/rich_change_',cases[c], '.png'), 
  #        dpi=600,
  #        width = 14, height = 16, scale=1.5, units = 'cm'
  # )
  # 
  # dir.create('figures/svg_rich_change', showWarnings = F)
  # # Save
  # ggsave(plot = final_plot, 
  #        filename = paste0('figures/svg_rich_change/', 
  #                          cases[c], '.svg'), 
  #        dpi=600,
  #        width = 14, height = 16, scale=1.2, units = 'cm'
  # )
}

rastOut <- rastOut %>% rast()
names(rastOut) <- cases

writeRaster(disagg(rastOut, fact=8, method='bilinear') %>% mask(ca), 'figures/plants_all_future_rich_changes.tif')

rastRichFut <- rastRichFut %>% rast()
names(rastRichFut) <- cases

writeRaster(rastRichFut, 'figures/plants_all_future_rich_raw.tif')

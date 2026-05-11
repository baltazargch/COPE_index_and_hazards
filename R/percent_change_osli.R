library(tidyverse)
library(terra)
library(sf)

curr <- rast('outputs/models/CA_osli_current_ensemble.tif') %>% mean()/1000

csv_mods <- read_csv('outputs/models/projections/projection_index.csv')

# Target group (must include presences & background class column)
tg_os <- readr::read_csv("inputs/records/pres_abs_oslig_target_group.csv",
                         show_col_types = FALSE)

# Presences to sf (EPSG:4326 assumed in CSV)
st_occs <- tg_os %>% filter(lon < -100) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)%>%
  filter(class != "background")


th <- terra::extract(curr, st_occs) %>% na.omit() %>% pull(2)
th <- quantile(th, 0.10) %>% round(3)

gr_mods <- csv_mods %>% 
  mutate(group = str_c(ssp,'_', period)) 

gr_mods <- gr_mods %>% split(.,.$group)

means_r <- map(gr_mods, ~(rast(.x$path)/1000) %>% mean()) %>% rast() %>% 
  project(curr) %>% crop(curr, mask=T)

plot(means_r - curr)

# writeRaster(means_r - curr, 'figures/osli_percent_changes.tif')

map(1:9, \(.x) 
    {
      ggplot(tibble(vals = values(means_r[[.x]]))) + 
        geom_histogram(aes(x=vals, fill=vals>th), binwidth = 0.01) + 
        ggsci::scale_fill_uchicago(alpha = 0.6) +
        # xlim(c(0,1))+
        geom_vline(xintercept = th, linetype=2, linewidth=0.55, colour='gray30') +
        labs(title = names(means_r)[.x], 
             x = 'Suitability', y='Frequency')+
        theme_light(base_size = 10) 
      
      ggsave(paste0('figures/suit_osli_percent/', names(means_r)[.x], '_hist.png'))
})

#extract table of areas
all_r <- c(curr, means_r)

names(all_r)[1] <- 'current'

df_areas <- map_dfr(1:nlyr(all_r), 
                    \(x){
                      # x <- all_r[[1]]
                      zns <- all_r[[x]] > th
                      
                      plot(zns)
                      
                      areas <- expanse(zns, byValue = T, unit='km')
                      areas %>% 
                        mutate(percent = round((area / sum(.data$area) *100), 2)) %>% 
                        mutate(layer = names(all_r[[x]]))
                    })

df_pub <- df_areas %>% 
  # split "layer" into scenario + period
  separate(layer, into = c("scenario_raw", "period"), sep = "_", fill = "right") %>% 
  mutate(
    # nice scenario labels
    scenario = recode(
      scenario_raw,
      "current" = "Current",
      "ssp245" = "SSP2-4.5",
      "ssp370" = "SSP3-7.0",
      "ssp585" = "SSP5-8.5"
    ),
    # for current, period is just "Current"
    period   = if_else(is.na(period), "1950-2014", period),
    # recode suitability
    suitability = if_else(value == 1, "Suitable", "Unsuitable"),
    # clearer column names
    `Area (km²)` = area,
    `Area (%)`   = percent
  ) %>% 
  select(scenario, period, suitability, `Area (km²)`, `Area (%)`) %>% 
  mutate(
    scenario = factor(scenario, levels = c("Current", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")),
    period   = factor(period, levels = c("1950-2014", "2015-2044", "2045-2074", "2075-2100")),
    suitability = factor(suitability, levels = c("Suitable", "Unsuitable"))
  ) %>% 
  arrange(scenario, period) 
df_pub_wide <- df_pub %>% 
  pivot_wider(
    names_from = suitability,
    values_from = c(`Area (km²)`, `Area (%)`),
    names_glue = "{suitability} { .value }"
  ) %>% 
  arrange(scenario, period)

# World basemap via rnaturalearth (land polygons)
wrld_org <- rnaturalearth::ne_states(country = 'united states of america', 
                                     returnclass = "sf")

ca <- wrld_org %>% filter(name_en == 'California')

#base layer to extend when plants are not in whole California
r <-  rast("outputs/models/NA_osli_current_ensemble.tif") %>% 
  crop(ca, mask=T) / 1000
wrld <- wrld_org %>% 
  st_crop(extend(ext(r), 0.5))

wrldp <- st_transform(wrld, crs=4269)
walk(1:nlyr(means_r), \(x){
  # x = 1
  # Set all to desired CSR
  difrp <- disagg(100*(means_r[[x]] - curr), fact=8, method='bilinear') %>% 
    mask(ca) %>% project(crs('epsg:4269')) 
  
  # make the color range 
  # Keep your actual range (no symmetry), and label ticks as percentages
  rng   <- as.numeric(terra::global(difrp, fun = "range", na.rm = TRUE))   # c(min, max)
  brks  <- pretty(rng, n = 6)
  
  # --- (B) Main map (your plot, lightly tidied) -------------------------------
  p_map <- ggplot() +
    geom_sf(data = wrldp, colour = "gray50", fill = "gray88", show.legend = FALSE) +
    tidyterra::geom_spatraster(data = difrp, 
                               interpolate = T, alpha = 0.9, na.rm = TRUE) +
    scale_fill_gradient2(
      low = 'darkred', 
      mid = "#F2F2F2", 
      high = 'darkgreen',
      midpoint = 0,
      # limits = c(-100, 100),
      labels = scales::label_number(suffix = "%"), name = "Suitability change",
      oob = scales::squish, na.value = NA
    )+
    coord_sf(expand = F, xlim = c(-124.5, -114), ylim = c(32.5, 42.2), crs = 4269) +
    ggspatial::annotation_scale(text_cex=1,
                                location = "bl", width_hint = 0.25, style = "ticks",
                                line_col = "gray20", text_col = "gray20"
    ) +
    labs(x = NULL, y = NULL, 
         title = names(means_r)[x]
    ) +
    theme_minimal(base_size = 16) +
    theme(
      text = element_text(family = "sans"),
      panel.background = element_rect(fill = "#D9E8F5"),
      legend.position = "right",
      panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.5), 
    )
  
  # p_map
  # --- (C) Inset barplot ------------------------------------------------------
  # area-weighted bins (km²)
  cell_km2 <- terra::mask(terra::cellSize(difrp, unit="km"), difrp)
  df <- tibble(change = as.numeric(terra::values(difrp)),
               area = as.numeric(terra::values(cell_km2))) |>
    filter(is.finite(change), is.finite(area))
  
  breaks <- c(-Inf, -30, -5, 5, 30, Inf)
  
  labs <- c("≤ -30%",
            "-30 to -5%",
            "-5 to 5%",
            "5 to 30%",
            "≥ 30%")
  
  bin_df <- df |>
    mutate(bin = cut(change, breaks = breaks, labels = labs, include.lowest = TRUE)) |>
    count(bin, wt = area, name = "area_km2") |>
    tidyr::complete(bin, fill = list(area_km2 = 0)) |>
    mutate(pct = 100 * area_km2 / sum(area_km2),
           bin = factor(bin, levels = labs))
  
  p_inset <- ggplot(bin_df, aes(x = bin, y = pct, fill = bin)) +
    geom_col(width = 0.8, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.0f%%", pct)), vjust = -0.25, size = 3) +
    scale_fill_manual(values = c("#B8574E","#E6B0A4","#F2F2F2","#A7C7B7","#228833")) +
    scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
    labs(x = NULL, y = "Area (%)") +
    theme_light(base_size = 10) +
    theme(
      plot.margin = margin(4, 4, 4, 4),
      axis.title.y = element_text(margin = margin(r = 2)),
      # axis.text.x = element_text(angle = 0, vjust = 1),
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(),
      plot.background = element_rect(color='gray30', fill = 'white'),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )
  
  
  # --- (D) Compose with patchwork (top-right inset) ---------------------------
  # Adjust the numbers [0..1] to fine-tune placement if needed.
  final_plot <- p_map + patchwork::inset_element(
    p_inset,
    left   = 0.63,  # x0
    bottom = 0.62,  # y0
    right  = 0.99,  # x1
    top    = 0.98,  # y1
    align_to = "plot"
  ) 
  
  # final_plot
  
  # Save
  ggsave(plot = final_plot, 
         filename = paste0('figures/suit_osli_percent/png/', names(means_r)[x], '.png'), 
         dpi=600,
         width = 14, height = 16, scale=1.5, units = 'cm'
  )
  
  # Save
  ggsave(plot = final_plot, 
         filename = paste0('figures/suit_osli_percent/svg/', names(means_r)[x], '.svg'), 
         dpi=600,
         width = 14, height = 16, scale=1.2, units = 'cm'
  )
})

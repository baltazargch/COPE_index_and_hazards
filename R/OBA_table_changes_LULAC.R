# Load libraries ----------------------------------------------------------
library(terra)
library(sf)
library(tidyverse)
library(scales)

# Dictionary for LULAC model ----------------------------------------------
lulac <- tibble::tribble(
  ~value,                      ~landuse,
  1L,                          "Water",
  2L,                      "Developed",
  3L, "Transportation/Other Developed",
  4L,              "Mining (not Used)",
  5L,                         "Barren",
  6L,                         "Forest",
  7L,                      "Grassland",
  8L,             "Annual Agriculture",
  9L,                       "Wetlands",
  10L,                     "Shrublands",
  11L,                        "SnowIce",
  12L,          "Perennial Agriculture"
)

# Load and clean data on land use -----------------------------------------
lu_hist <- list.files('/mnt/4TB/GIS/Rasters/LULAC California BAU/Historical/StateClass_Historical_1970-2001/State Class/', 
                      '.tif', full.names = T)

lu_low <- list.files('/mnt/4TB/GIS/Rasters/LULAC California BAU/BAU_Low/StateClass_Projection_BAU+Low_2001-2101/State Class Projection PopLow - 2001-2100/', 
                     '.tif', full.names = T)

lu_high <- list.files('/mnt/4TB/GIS/Rasters/LULAC California BAU/BAU_High/StateClass_Projection_BAU+High_2001-2101/State Class Projection PopHigh - 2001-2100/', 
                      '.tif', full.names = T)


#  UDF for repetitive actions  --------------------------------------------
# Extracting year from basename of files 
extract_idx <- function(x, low=NULL, high=NULL){
  # x <- lu_low
  x <- x %>% 
    basename() %>% 
    str_extract('(?<=Ts)[0-9]{4}(?=-SClass\\.tif)')
  
  between(x, low, high)
}

# Calculate modal and name raster
load_raster <- function(x, name = ''){
  x <- x %>% 
    rast() %>% 
    modal() %>% 
    as.numeric()
  set.names(x, name)
  x
}

# Load and name suitability for Osmia and divide if not cloglog output
load_models <- function(x, ref, name = ''){
  x <- x %>% rast() %>% 
    project(ref) %>%
    crop(ref, mask=T) %>% 
    median()
  set.names(x, name)
  if(max(minmax(x)) > 1) x / 1000 else x 
}

# Historical data is from 1950 - 2001 so it is capped
lu_hist_rast <- load_raster(lu_hist[extract_idx(lu_hist, 1990, 2002)], 
                            'lu_hist')

# Load end decade modal for each population trend in BAU ------------------
# BAU low
lu_low_2044_rast <- load_raster(lu_low[extract_idx(lu_low, 2034, 2044)], 
                                'lu_low_2044')

lu_low_2074_rast <- load_raster(lu_low[extract_idx(lu_low, 2064, 2074)], 
                                'lu_low_2074')

lu_low_2100_rast <- load_raster(lu_low[extract_idx(lu_low, 2090, 2100)], 
                                'lu_low_2100')

# BAU high
lu_high_2044_rast <- load_raster(lu_high[extract_idx(lu_high, 2034, 2044)], 
                                 'lu_high_2044')

lu_high_2074_rast <- load_raster(lu_high[extract_idx(lu_high, 2064, 2074)], 
                                 'lu_high_2074')

lu_high_2100_rast <- load_raster(lu_high[extract_idx(lu_high, 2090, 2100)], 
                                 'lu_high_2100')

# Load Osmia projections  -------------------------------------------------
osli_hist_rast <- load_models('outputs/models/NA_osli_current_ensemble.tif',
                              lu_hist_rast, 
                              'current_suitabilty')  

osli_2044_bau <- load_models('outputs/models/projections/means/ssp585_2015-2044_EMmean.tif',
                             lu_hist_rast, 
                             '2044_suitabilty')

osli_2074_bau <- load_models('outputs/models/projections/means/ssp585_2045-2074_EMmean.tif',
                             lu_hist_rast, 
                             '2074_suitabilty')

osli_2100_bau <- load_models('outputs/models/projections/means/ssp585_2075-2100_EMmean.tif',
                             lu_hist_rast, 
                             '2100_suitabilty')

# Calculate areas ---------------------------------------------------------
## UDF area ----------------------------------------------------------------
calculate_period_area <- function(lu_base, 
                                  lu_period,
                                  suit_base,
                                  suit_period, 
                                  lulac, 
                                  period_label){
  # lu_base = lu_hist_rast
  # lu_period = lu_high_2100_rast
  # suit_base = osli_hist_rast
  # suit_period = osli_2100_bau
  # period_label = '2100'
  
  msk <- !is.na(lu_base) & !is.na(lu_period) & !is.na(suit_base) & !is.na(suit_period)
  lc_b <- mask(lu_base,  msk)
  lc_f <- mask(lu_period, msk)
  Sb   <- mask(suit_base,   msk)
  Sf   <- mask(suit_period, msk)
  
  # cell area (km2)
  a_km2 <- cellSize(lc_b, unit = "km") %>% mask(lc_b)
  
  # transitions
  loss <- (Sf - Sb) < 0    # loss of suitability
  gain <- (Sf - Sb) > 0    # gain of suitability
  
  # losses/gains by COVER (stratify by FUTURE LULC for the period/trend)
  loss_by_cover <- zonal( (loss * a_km2), lc_f, fun = "sum", na.rm = TRUE) |>
    as_tibble() |>
    rename_with(~c("value","loss_km2"), 1:2)
  
  gain_by_cover <- zonal( (gain * a_km2), lc_f, fun = "sum", na.rm = TRUE) |>
    as_tibble() |>
    rename_with(~c("value","gain_km2"), 1:2)
  
  df <-  full_join(loss_by_cover, gain_by_cover, by = "value") %>% 
    left_join(lulac, by = "value") %>% 
    mutate(period_label=period_label) %>% 
    select(landuse, value, loss_km2, gain_km2, period_label)
  df
}

lu_by_trend <- list(
  low  = list(`2044` = lu_low_2044_rast,
              `2074` = lu_low_2074_rast,
              `2100` = lu_low_2100_rast),
  high = list(`2044` = lu_high_2044_rast,
              `2074` = lu_high_2074_rast,
              `2100` = lu_high_2100_rast)
)

suit_by_period <- list(
  `2044` = osli_2044_bau,
  `2074` = osli_2074_bau,
  `2100` = osli_2100_bau
)

cases <- tidyr::expand_grid(
  period = c("2044","2074","2100"),
  bau    = c("low","high")
) %>%
  mutate(case = str_c(bau, period, sep = "_"))

tb_list <- purrr::pmap(cases, \(period, bau, case) {

  lu_period    <- lu_by_trend[[bau]][[period]]
  suit_period  <- suit_by_period[[period]]
  
  # error managing
  stopifnot(!is.null(lu_period), !is.null(suit_period))
  
  calculate_period_area(
    lu_base      = lu_hist_rast,
    lu_period    = lu_period,
    suit_base    = osli_hist_rast,
    suit_period  = suit_period,
    lulac        = lulac,
    period_label = case
  )
})


# bind into one table
table_all <- dplyr::bind_rows(tb_list)


table_wide <- table_all %>%
  mutate(period = str_extract(period_label, '[0-9]{4}')) %>% 
  mutate(population_trend = str_remove(period_label, '_.*')) %>% 
  arrange(landuse, population_trend, period) %>%
  pivot_wider(
    id_cols = c(landuse, population_trend),
    names_from = period,
    values_from = c(loss_km2, gain_km2)
  ) %>%
  arrange(landuse, population_trend) %>%
  mutate(landuse = case_when(
    str_detect(landuse, 'Agriculture') ~ "Agriculture", 
    str_detect(landuse, 'Developed') ~ "Developed",
    .default = landuse 
  )) %>% 
  group_by(landuse, population_trend) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)), .groups = "drop")

write_csv(table_wide, 'outputs/csv/area_per_poptrend_period_landuse.csv')


# Plot it! ----------------------------------------------------------------
# glimpse(table_wide)

df_long <- table_wide  %>% 
  filter(
    !landuse  %in% c('Water', 'SnowIce', 'Barren', 'Wetlands')
  ) %>% 
  pivot_longer(
    cols = matches("^(loss|gain)_km2_\\d+$"),
    names_to = c("type", ".value", "period"),
    names_pattern = "^(loss|gain)_(km2)_(\\d+)$"
  ) %>%
  # km2 now holds the numeric area; make losses negative
  mutate(
    km2_signed = if_else(type == "loss", -km2, km2),
    period = factor(period, levels = c("2044","2074","2100")),
    population_trend = factor(population_trend, levels = c("low","high")),
    type = factor(type, levels = c("loss","gain"), labels = c("Area loss","Area gain"))
  )

# (optional) order land uses by overall magnitude so facets read nicely
lu_order <- df_long %>%
  group_by(landuse) %>%
  summarise(mag = sum(abs(km2_signed), na.rm = TRUE), .groups = "drop") %>%
  arrange(mag) %>% pull(landuse)

df_long <- df_long %>% mutate(landuse = factor(landuse, levels = lu_order)) %>% 
  filter(!landuse  %in% c('Barren', '')) %>% 
  filter(population_trend == 'high')

# --- Step 1. Compute totals for each land use ----------------------------
cover_totals <- table_wide %>%
  group_by(landuse) %>%
  summarise(total_area = sum(
    loss_km2_2044 + loss_km2_2074 + loss_km2_2100 +
      gain_km2_2044 + gain_km2_2074 + gain_km2_2100,
    na.rm = TRUE
  ))

# --- Step 2. Reshape to long and join totals -----------------------------
df_long <- table_wide %>%
  filter(
    !landuse  %in% c('Water', 'SnowIce', 'Barren', 'Wetlands')
  ) %>% 
  pivot_longer(
    cols = starts_with("loss_km2_") | starts_with("gain_km2_"),
    names_to = c("type","period"),
    names_pattern = "(loss|gain)_km2_(\\d+)",
    values_to = "km2"
  ) %>%
  mutate(
    km2_signed = ifelse(type == "loss", -km2, km2)
  ) %>%
  left_join(cover_totals, by = "landuse") %>%
  mutate(
    rel_change = km2_signed / total_area   # relative to cover size
  )

# 2) Plot: red (loss) below zero, dark green (gain) above; facet by trend x period
df_long %>% 
  filter(population_trend == 'high') %>% 
  ggplot( aes(x = rel_change, y = landuse, fill = type)) +
  geom_col(position = "identity", alpha=0.85) +
  geom_text(aes(label = scales::percent(rel_change, accuracy = 0.1)),
            hjust = 1,
            size = 2.8) +
  facet_grid(population_trend ~ period) +
  scale_fill_manual(values = c("loss" = "#B22222", "gain" = "#006400")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Relative area change by land use",
    subtitle = "Relative to total cover area; gains positive, losses negative",
    x = "Relative change (% of land cover area)",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", 
        panel.border = element_rect(colour = 'gray20'),
        axis.text.x = element_text(angle = 45, hjust = 1)
        )

plot(lu_high_2044_rast, type='classes')
plot(lu_high_2100_rast, type='classes')

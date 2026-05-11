library(terra)
library(sf)
library(tidyverse)
library(tidyterra) 
library(scales)

lu_hist <- list.files('/mnt/4TB/GIS/Rasters/LULAC California BAU/Historical/StateClass_Historical_1970-2001/State Class/', 
                      '.tif', full.names = T)

lu_bau <- list.files('/mnt/4TB/GIS/Rasters/LULAC California BAU/BAU/StateClass_Projection_BAU_2001-2101/State Class/', 
                     '.tif', full.names = T)

lu_bau_high <- list.files('/mnt/4TB/GIS/Rasters/LULAC California BAU/BAU_High/StateClass_Projection_BAU+High_2001-2101/State Class Projection PopHigh - 2001-2100/', 
                          '.tif', full.names = T)


idx  <- lu_hist %>% 
  basename() %>% 
  str_extract('(?<=Ts)[0-9]{4}(?=-SClass\\.tif)') %>%
  as.numeric()

lu_hist_rast <- rast(lu_hist[which(idx > 1979)])

idx  <- lu_bau %>% 
  basename() %>% 
  str_extract('(?<=Ts)[0-9]{4}(?=-SClass\\.tif)') %>%
  as.numeric()

lu_bau_rast <- rast(lu_bau[between(idx, 2014, 2014)])

idx  <- lu_bau_high %>% 
  basename() %>% 
  str_extract('(?<=Ts)[0-9]{4}(?=-SClass\\.tif)') %>%
  as.numeric()

lu_bau_high_rast <- rast(lu_bau_high[idx > 2099])


lu_hist_rast <- modal(lu_bau_rast) %>% as.numeric()
names(lu_hist_rast) <- 'historical'

lu_bau_rast <- modal(lu_bau_high_rast) %>% as.numeric()
names(lu_bau_rast) <- 'projected'


osli_suit_curr <-osli_suit_curr <- rast('outputs/models/NA_osli_current_ensemble.tif')%>% 
  project(lu_hist_rast) %>%
  crop(lu_hist_rast, mask=T) %>% 
  median()/1000
names(osli_suit_curr) <- 'suitability_current'

osli_suit_fut <- rast("outputs/models/projections/means/ssp585_2075-2100_EMmean.tif") %>% 
  project(lu_hist_rast) %>%
  crop(lu_hist_rast, mask=T) %>% 
  median()/1000
names(osli_suit_fut) <- 'suitability_future'

db_curr <- cbind(
  as.data.frame(lu_hist_rast, na.rm=F),
  as.data.frame(lu_bau_rast, na.rm=F),
  as.data.frame(osli_suit_curr, na.rm=F),
  as.data.frame(osli_suit_fut, na.rm=F)
) %>% as_tibble() %>% filter(if_all(everything(), ~ !is.na(.)))

lulac <- tibble::tribble(
  ~Value,                      ~LandUse,
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

plot(lu_hist_rast == 11)

# plot(lu_hist_rast == 7)

esm_osli <- readRDS('outputs/models/osmia_bm_ensemble_model.rds')
evals <- get_evaluations(esm_osli)
tss_thresh_ens <- evals %>%
  filter(metric.eval == "TSS") %>%
  pull(cutoff) %>% median()/1000


df_long <- db_curr %>%
  # Map codes → land-use names
  mutate(
    LandUse_hist = case_when(
      historical %in% c(2,3)         ~ "Developed",
      historical %in% c(6) ~ "Forest",
      historical %in% c(7) ~ "Grassland",
      historical %in% c(9) ~ "Wetlands",
      historical %in% c(10) ~ "Shrublands",
      historical %in% c(8,12)        ~ "Agriculture",
      .default = NA
    ),
    LandUse_proj = case_when(
      projected %in% c(2,3)         ~ "Developed",
      projected %in% c(6) ~ "Forest",
      projected %in% c(7) ~ "Grassland",
      projected %in% c(9) ~ "Wetlands",
      projected %in% c(10) ~ "Shrublands",
      projected %in% c(8,12)        ~ "Agriculture",
      .default = NA
    )
  ) %>%
  # pivot longer on suitability columns
  pivot_longer(
    cols = c(suitability_current, suitability_future),
    names_to = "Period",
    values_to = "Suitability"
  ) %>%
  mutate(
    Period = recode(Period,
                    "suitability_current" = "Historical",
                    "suitability_future"  = "Projected"
    ),
    # pick correct land-use according to Period
    LandUse = if_else(Period == "Historical", LandUse_hist, LandUse_proj),
    Period  = factor(Period, levels = c("Historical", "Projected")),
    LandUse = factor(LandUse, levels = names(lu_cols))
  ) %>%
  filter(LandUse_hist != 'Barren/Other') %>% 
  filter(!is.na(LandUse))

# land use palette
lu_cols <- c(
  "Developed"     = "#e34a33",
  "Forest"       = "darkgreen",
  "Grassland"       = "#31a984",
  "Wetlands"       = "#3182bd",
  "Shrublands"       = "#8c8c39",
  "Agriculture"   = "#fdae61"
)

# ---- Plot: same hue, darker for Projected ----
# ---- Bin suitability and count ----
binwidth <- 0.06   # adjust as needed

s_min <- floor(min(df_long$Suitability) / binwidth) * binwidth
s_max <- ceiling(max(df_long$Suitability) / binwidth) * binwidth
breaks <- seq(s_min, s_max, by = binwidth)

df_binned <- df_long %>%
  mutate(
    bin = cut(Suitability, breaks = breaks, include.lowest = TRUE, right = FALSE),
    bin_mid = (as.numeric(bin) - 0.5) * binwidth + s_min
  ) %>%
  count(LandUse, Period, bin_mid, name = "n") %>%
  group_by(LandUse, Period) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(df_binned, aes(x = bin_mid, y = prop,
                      fill = LandUse, alpha = Period)) +
  geom_col(position = position_dodge(width = binwidth * 0.9),
           color = "white", width = binwidth * 0.9) +
  facet_wrap(LandUse~.,ncol=2, scales='free') +
  scale_fill_manual(values = lu_cols) +
  scale_alpha_manual(values = c(Historical = 0.45, Projected = 0.9)) +
  scale_x_continuous(breaks = pretty_breaks(7)) +
  labs(
    x = "Suitability",
    y = "Proportion within Period",
    fill = "Land Use",
    alpha = "Period",
    title = "Suitability Distributions by Land Use\nHistorical vs Projected"
  ) +
  envalysis::theme_publish(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) + 
  geom_vline(xintercept = tss_thresh_ens, linetype = 2) 

ggsave("figures/lu_svgs/distribution_suitability_lu.svg", 
       width = 6, height = 4, dpi = 300, scale = 2)

df_long  %>%
  mutate(
    Period_group = str_c(Period, LandUse)
  ) %>% 
  ggplot(aes(x = LandUse, y = Suitability,
             fill = LandUse, group = Period_group, alpha=Period)) +
  geom_boxplot(outliers = F) +
  # geom_violin() +
  scale_fill_manual(values = lu_cols) +
  scale_alpha_manual(values = c(Historical = 0.9, Projected = 0.45)) +
  labs(
    x = "Land Use",
    y = "Suitability",
    alpha = "Period",
    fill = "Land Use",
    title = "Suitability by Land Use: Historical (lighter) vs Projected (darker)"
  ) +
  envalysis::theme_publish(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )
ggsave("figures/lu_svgs/boxplot_suitability_lu.svg", 
       width = 6, height = 4, dpi = 300, scale = 2)
# Map land-use names to the underlying codes in the raster
lu_codes <- list(
  "Developed"     = c(2, 3),
  "Water"         = c(1),
  "Barren/Other"  = c(4, 5),
  "Natural"       = c(6, 7, 9, 10, 11),
  "Agriculture"   = c(8, 12)
)

# Matching, intuitive colors
lu_cols <- c(
  "Developed"     = "#e34a33",  # red/orange
  "Water"         = "#3182bd",  # blue
  "Barren/Other"  = "#969696",  # gray
  "Natural"       = "#31a354",  # green
  "Agriculture"   = "#fdae61"   # yellow/orange
)

dir.create("figures/lu_svgs", showWarnings = FALSE)


iwalk(lu_codes, function(codes, nm) {
  # codes = lu_codes[[1]]
  # nm = names(lu_codes)[1]
  # Binary raster: 1 = class cells, 0 = others
  r_bin <- ifel(lu_hist_rast %in% codes, 1, 0) %>% mask(lu_hist_rast)
  names(r_bin) <- "class_bin"
  
  
  # Build plot
  p <- ggplot() +
    geom_spatraster(data =as.factor(r_bin)) +
    scale_fill_manual(na.value = NA,
                      values = c("gray90", lu_cols[[nm]]),
                      guide = "none"
    ) +
    coord_sf() +
    theme_void()
  
  # Save one SVG per land use
  out <- file.path("figures/lu_svgs", 
                   tolower(paste0("current_", str_replace_all(nm, "/", "-"))))
  ggsave(paste0(out, ".svg"), p, width = 6, height = 4, dpi = 300)
})
iwalk(lu_codes, function(codes, nm) {
  # codes = lu_codes[[1]]
  # nm = names(lu_codes)[1]
  # Binary raster: 1 = class cells, 0 = others
  r_bin <- ifel(lu_bau_rast %in% codes, 1, 0) %>% mask(lu_hist_rast)
  names(r_bin) <- "class_bin"
  
  
  # Build plot
  p <- ggplot() +
    geom_spatraster(data =as.factor(r_bin)) +
    scale_fill_manual(na.value = NA,
                      values = c("gray90", lu_cols[[nm]]),
                      guide = "none"
    ) +
    coord_sf() +
    theme_void()
  
  # Save one SVG per land use
  out <- file.path("figures/lu_svgs", 
                   tolower(paste0("future_", str_replace_all(nm, "/", "-"))))
  ggsave(paste0(out, ".svg"), p, width = 6, height = 4, dpi = 300)
})

r_lu <- app(lu_hist_rast, fun = function(v) {
  ifelse(v %in% c(1),                      1,  # Water
         ifelse(v %in% c(2, 3),                   2,  # Developed
                ifelse(v %in% c(4, 5),                   3,  # Barren/Other
                       ifelse(v %in% c(6, 7, 9, 10, 11),        4,  # Natural
                              ifelse(v %in% c(8, 12),                  5,  # Agriculture
                                     NA)))))
})

# 2) Make it a categorical raster with labels in the right order
r_lu <- as.factor(r_lu)
lev <- levels(r_lu)[[1]]
colnames(lev) <- tolower(colnames(lev))

lev$label <- c("Water", "Developed", "Barren/Other", "Natural", "Agriculture")[match(lev$id, 1:5)]
levels(r_lu) <- as.data.frame(lev[,c(1,3)])
r_lu <- droplevels(r_lu) 
# 3) Your palette (same hues as before)
lu_cols <- c(
  "Developed"     = "#e34a33",
  "Water"         = "#3182bd",
  "Barren/Other"  = "#969696",
  "Natural"       = "#31a354",
  "Agriculture"   = "#fdae61"
)

# 4) Plot all classes in one map
p <- ggplot() +
  geom_spatraster(data = r_lu) +
  scale_fill_manual(
    values   = lu_cols,
    na.value = NA,
    drop     = TRUE,
    name     = "Land Use"
  ) +
  coord_sf() +
  theme_void()
p
ggsave("figures/lu_svgs/map_lu_current.svg", 
       width = 6, height = 4, dpi = 300, scale = 2)

r_lu <- app(lu_bau_rast, fun = function(v) {
  ifelse(v %in% c(1),                      1,  # Water
         ifelse(v %in% c(2, 3),                   2,  # Developed
                ifelse(v %in% c(4, 5),                   3,  # Barren/Other
                       ifelse(v %in% c(6, 7, 9, 10, 11),        4,  # Natural
                              ifelse(v %in% c(8, 12),                  5,  # Agriculture
                                     NA)))))
})

# 2) Make it a categorical raster with labels in the right order
r_lu <- as.factor(r_lu)
lev <- levels(r_lu)[[1]]
colnames(lev) <- tolower(colnames(lev))

lev$label <- c("Water", "Developed", "Barren/Other", "Natural", "Agriculture")[match(lev$id, 1:5)]
levels(r_lu) <- as.data.frame(lev[,c(1,3)])
r_lu <- droplevels(r_lu) 
# 3) Your palette (same hues as before)
lu_cols <- c(
  "Developed"     = "#e34a33",
  "Water"         = "#3182bd",
  "Barren/Other"  = "#969696",
  "Natural"       = "#31a354",
  "Agriculture"   = "#fdae61"
)

# 4) Plot all classes in one map
p <- ggplot() +
  geom_spatraster(data = r_lu) +
  scale_fill_manual(
    values   = lu_cols,
    na.value = NA,
    drop     = TRUE,
    name     = "Land Use"
  ) +
  coord_sf() +
  theme_void()
p
ggsave("figures/lu_svgs/map_lu_bau.svg", 
       width = 6, height = 4, dpi = 300, scale = 2)



# your mapping and palette
lu_codes <- list(
  "Developed"     = c(2, 3),
  "Water"         = c(1),
  "Barren/Other"  = c(4, 5),
  "Natural"       = c(6, 7, 9, 10, 11),
  "Agriculture"   = c(8, 12)
)
lu_cols <- c(
  "Developed"     = "#e34a33",
  "Water"         = "#3182bd",
  "Barren/Other"  = "#969696",
  "Natural"       = "#31a354",
  "Agriculture"   = "#fdae61"
)

# 1) Make a tidy lookup table: code -> LandUse
lu_lookup <- stack(lu_codes) |>
  dplyr::rename(code = values, LandUse = ind)   # columns: code, LandUse

# 2) Build the long df and join the lookup
df_area_both <- db_curr %>%
  select(historical, projected) %>%
  pivot_longer(everything(), names_to = "Period", values_to = "code") %>%
  mutate(Period = recode(Period, historical = "Historical", projected = "Projected")) %>%
  left_join(lu_lookup, by = "code") %>%        # <-- no vapply; fast & clean
  filter(!is.na(LandUse)) %>%
  count(Period, LandUse, name = "n") %>%
  group_by(Period) %>%
  mutate(Percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    LandUse = factor(LandUse, levels = names(lu_cols)),
    Period  = factor(Period, levels = c("Historical","Projected"))
  )

# 3) Plot: dodged percent bars (Projected darker via alpha)
ggplot(df_area_both, aes(x = LandUse, y = Percent, fill = LandUse, alpha = Period)) +
  geom_col(position = position_dodge(width = 0.8), color = "white", width = 0.75) +
  scale_fill_manual(values = lu_cols) +
  scale_alpha_manual(values = c(Historical = 0.9, Projected = 0.45)) +
  labs(
    x = "Land Use", y = "Percent of Total Area (%)",
    title = "Relative Area by Land Use: Historical vs Projected",
    alpha = "Period", fill = "Land Use"
  ) +
  envalysis::theme_publish(base_size = 14)+ 
  geom_text(aes(label = sprintf("%.1f%%", Percent)),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3)+
  
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave("figures/lu_svgs/relative_lu_area.svg", 
       width = 6, height = 4, dpi = 300, scale = 2)



######### pedning ######

library(terra)
library(sf)
library(tidyverse)

# World basemap via rnaturalearth (land polygons)
wrld_org <- rnaturalearth::ne_states(country = 'united states of america', 
                                     returnclass = "sf")

ca <- wrld_org %>% filter(name_en == 'California') %>% 
  st_transform(3310) #NAD83

osli_ssp245_t1 <- rast('outputs/models/projections/means/ssp245_2015-2044_EMmean.tif') %>% 
  project(ca) %>% 
  crop(ca, mask=T)/1000

all_fls <- list.files('outputs/final_maxent/rasters/', 'p10.tif$', 
                      full.names = T)

current_richness <- all_fls[ !str_detect(all_fls, 'ssp') ]

c_rich <- map(current_richness, ~ rast(.x) %>% 
                project(ca) %>%
                crop(ca, mask=T) %>% 
                resample(osli_ssp245_t1, method='near')) 

curr_resource_rich <- c_rich %>% 
  rast() %>% 
  sum(na.rm = T) %>% 
  mask(osli_ssp245_t1)

future_richness <- all_fls[ str_detect(all_fls, '2015-2044_ssp245') ]

f_rich <- map(future_richness, ~ rast(.x) %>% 
                project(ca) %>%
                crop(ca, mask=T) %>% 
                resample(osli_ssp245_t1, method='near')) 

future_resource_rich <- f_rich %>% 
  rast() %>% 
  sum(na.rm = T) %>% 
  mask(osli_ssp245_t1)

dR_log  <- log((future_resource_rich + 1) / (curr_resource_rich + 1))
dR_tanh <- tanh(dR_log)  # smoothly bounded [-1,1]

w     <- (dR_tanh + 1) / 2 

# --- Power-tilt parameters ---
lambda <- 0.7  # 0 < lambda <= 0.5; higher = stronger boost/damp

# alpha(w) = 1 - lambda*(2w - 1)
alpha <- 1 - lambda * (2*w - 1)

# Composite index: C = S ^ alpha(w)
C <- osli_ssp245_t1 ^ alpha 

library(patchwork)

plot(alpha)

plot(C) 

plot(osli_ssp245_t1)


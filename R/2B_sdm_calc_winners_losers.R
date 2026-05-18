library(tidyverse)
library(scales)
library(terra)
library(sf)
library(rnaturalearth)

# -----------------------------
# 1) Read TSS cutoffs
# -----------------------------
evals <- read_csv("outputs/csv/plants_ensemble_eval.csv") %>% 
  filter(metric.eval == "TSS") %>% 
  mutate(
    species = str_remove(full.name, "_.*") %>% 
      str_replace_all("\\.", "_"),
    .before = 1
  ) %>% 
  distinct(species, cutoff)

# ---- Load plants occs and list ----
plant_list <- read_csv('outputs/records/all_species_records_and_native.csv')
natives <- plant_list %>% filter(`CA Native` == 'yes')

fls_occs <- list.files('outputs/records/plants', '.csv', full.names = T)

natives$path <- sapply(natives$species, \(x) fls_occs[grep(str_replace(x, ' ', '_'), fls_occs)], 
                       simplify = T) %>% unname()
natives$sp <- natives$species %>% str_to_lower() %>% str_replace(' ', '_')

# -----------------------------
# 2) Read projection table
# -----------------------------
csv_mods <- read_csv("outputs/models/plants_projection_table.csv") %>% 
  mutate(
    dir_path = str_remove(path, basename(path)), 
    binary_file = path %>% basename() %>% sub("^([^_]+)_(.+\\.tif)$", "\\1_TSSbinay_\\2", .) %>% str_c(dir_path, .), 
    exists = file.exists(binary_file),
    current_file = paste0("outputs/models/", species, "_NA_cur_ensemble.tif")
  ) %>% 
  left_join(evals, by = "species") %>% 
  left_join(natives %>% select(species= sp, occs_path=path), 
            by = 'species')

# -----------------------------
# 3) California polygon
# -----------------------------
wrld_org <- rnaturalearth::ne_states(
  country = "united states of america",
  returnclass = "sf"
)

ca <- wrld_org %>% 
  filter(name_en == "California") %>% 
  st_make_valid()

# -----------------------------
# 4) Function to calculate suitable area
# -----------------------------
calc_binary_area <- function(r, cutoff, polygon_sf = ca) {
  
  polygon_v <- polygon_sf %>% 
    st_transform(crs(r)) %>% 
    vect()
  
  r_ca <- r %>% 
    crop(polygon_v) %>% 
    mask(polygon_v)
  
  suitable <- r_ca >= cutoff
  
  area <- expanse(suitable, unit = "km", byValue = TRUE) %>% 
    as_tibble()
  
  if (!any(area$value == 1)) {
    return(0)
  }
  
  area %>% 
    filter(value == 1) %>% 
    pull(area)
}

calc_baseline_info <- function(species, baseline_path, occs_path, polygon_sf = ca) {
  
  base_r <- rast(baseline_path)
  
  base_r <- if (nlyr(base_r) > 1) {
    mean(base_r, na.rm = TRUE)
  } else {
    base_r
  }
  
  occs <- read_csv(occs_path, show_col_types = FALSE)
  
  vals <- terra::extract(
    base_r,
    occs[, c("decimalLongitude", "decimalLatitude")]
  )
  
  cutoff <- vals %>% 
    select(-ID) %>% 
    pull(1) %>% 
    quantile(probs = 0.10, na.rm = TRUE) %>% 
    as.numeric()
  
  current_area_km2 <- calc_binary_area(
    r = base_r,
    cutoff = cutoff,
    polygon_sf = polygon_sf
  )
  
  tibble(
    species = species,
    cutoff = cutoff,
    current_area_km2 = current_area_km2
  )
}

baseline_info <- csv_mods %>% 
  distinct(species, current_file, occs_path) %>% 
  mutate(occs_path = ifelse(is.na(occs_path), "outputs/records/plants/Arctostaphylos_uva-ursi.csv", 
                          occs_path)) %>% 
  mutate(
    baseline_info = pmap(
      list(species, current_file, occs_path),
      \(species, current_file, occs_path) {
        calc_baseline_info(
          species = species,
          baseline_path = current_file,
          occs_path = occs_path,
          polygon_sf = ca
        )
      }
    )
  ) %>% 
  select(baseline_info) %>% 
  unnest(baseline_info)

baseline_info

calc_future_area <- function(future_path, cutoff, polygon_sf = ca) {
  
  fut_r <- rast(future_path)
  
  fut_r <- if (nlyr(fut_r) > 1) {
    mean(fut_r, na.rm = TRUE)
  } else {
    fut_r
  }
  
  calc_binary_area(
    r = fut_r,
    cutoff = cutoff,
    polygon_sf = polygon_sf
  )
}

areas_models <- csv_mods %>% 
  select(!cutoff) %>% 
  left_join(baseline_info, by = "species") %>% 
  mutate(
    future_area_km2 = map2_dbl(
      path,
      cutoff,
      \(path, cutoff) {
        calc_future_area(
          future_path = path,
          cutoff = cutoff,
          polygon_sf = ca
        )
      }
    ),
    area_change_km2 = future_area_km2 - current_area_km2,
    area_change_pct = if_else(
      current_area_km2 > 0,
      100 * area_change_km2 / current_area_km2,
      NA_real_
    )
  )


write_csv(
  areas_models,
  "outputs/plants_projected_area_change_10p.csv"
)


area_change_plot <- areas_models %>%
  filter(species != 'plagiobothrys_lithocaryus') %>% 
  mutate(
    period = factor(
      period,
      levels = c("2015-2044", "2045-2074", "2075-2100")
    )
  ) %>% 
  group_by(species, ssp, period) %>% 
  summarise(area_change_pct = mean(area_change_pct)) %>% 
  ungroup()



period_summary <- area_change_plot %>% 
  group_by(period, ssp) %>% 
  summarise(
    mean_change = mean(area_change_pct, na.rm = TRUE),
    sd_change = sd(area_change_pct, na.rm = TRUE),
    n = sum(!is.na(area_change_pct)),
    se_change = sd_change / sqrt(n),
    ci_low = mean_change - 1.96 * se_change,
    ci_high = mean_change + 1.96 * se_change,
    .groups = "drop"
  )


bin_breaks <- c(-Inf,-75,-50,-25,-10,10,25,50,75, Inf)
bin_labels <- c("≤ -75%","-75% to -50%","-50% to -25%","-25% to -10%",
                "-10% to 10%","10% to 25%","25% to 50%","50% to 75%","≥ 75%")

prep_counts <- function(df, col){
  df %>%
    mutate(change_bin = cut(.data[[col]], breaks = bin_breaks,
                            labels = bin_labels, ordered_result = TRUE)) %>%
    mutate(change_bin = forcats::fct_expand(change_bin, bin_labels)) %>%
    count(change_bin, name = "n_species", .drop = FALSE) %>% 
    mutate(Case = ifelse(str_detect(change_bin, '-'), 'Loss', 'Gain')) %>% 
    mutate(Case = ifelse(str_detect(change_bin, '-10% to'), 'No Change', Case)) %>% 
    mutate(Case = factor(Case, levels = c('Loss', 'No Change', 'Gain')))
}

low_2100 <- areas_models %>%
  filter(species != 'plagiobothrys_lithocaryus') %>% 
  filter(ssp == 'ssp245', period == '2075-2100') %>% 
  select(species, gcm, period, area_change_pct) %>% 
  prep_counts('area_change_pct')

high_2100 <- areas_models %>%
  filter(species != 'plagiobothrys_lithocaryus') %>% 
  filter(ssp == 'ssp585', period == '2075-2100') %>% 
  select(species, gcm, period, area_change_pct) %>% 
  prep_counts('area_change_pct')


rbind(low_2100 %>% mutate(Scenario = 'Low emissions', alpha='a') ,
      high_2100 %>% mutate(Scenario = 'High emissions', alpha='b')) %>% 
  mutate(Scenario = factor(Scenario, levels = c('Low emissions', 'High emissions'))) %>% 
  filter(!is.na(Case)) %>% 
  ggplot(aes(x = change_bin, y = n_species, fill = Case, alpha=Scenario)) +
  geom_bar(position = "dodge", stat='identity', color='white') +
  envalysis::theme_publish(base_size = 12) +
  scale_fill_manual(values=c('#B2182B', '#E6C85C', '#4DAF4A'))+
  guides(fill = "none")+
  scale_alpha_manual(values=c(0.5, 0.95))+
  labs(x = "% Change bin", y = "Number of species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = c(0.98, 0.98),      # top-right inside panel
        legend.justification = c(1, 1))

ggsave(filename = 'figures/win_los.png', 
       dpi=600,
       width = 9, height = 8, scale=1.5, units = 'cm')

area_change_plot_clean <- area_change_plot %>% 
  filter(species != 'plagiobothrys_lithocaryus') %>% 
  mutate(
    change_class = case_when(
      area_change_pct < -10 ~ "Loss",
      area_change_pct >  10 ~ "Gain",
      TRUE ~ "No change"
    ),
    change_class = factor(
      change_class,
      levels = c("Loss", "No change", "Gain")
    ),
    period = factor(
      period,
      levels = c("2015-2044", "2045-2074", "2075-2100")
    )
  ) %>% 
  filter(
    !is.na(area_change_pct)
  )

ggplot(area_change_plot_clean, aes(x = period, y = area_change_pct)) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_hline(
    yintercept = c(-10, 10),
    linewidth = 0.25,
    linetype = "dotted",
    colour = "grey60"
  ) +
  geom_violin(
    aes(group = period),
    width = 0.65,
    fill = "grey90",
    colour = "grey55",
    linewidth = 0.3,
    alpha = 0.45,
    trim = FALSE
  ) +
  geom_jitter(
    aes(shape = ssp, colour = change_class),
    width = 0.12,
    height = 0,
    alpha = 0.65,
    size = 1.9,
    stroke = 0.6
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 95,
    size = 8,
    linewidth = 0.9,
    colour = "black"
  ) +
  scale_y_continuous(
    labels = label_percent(scale = 1)
  ) +
  scale_colour_manual(
    values = c(
      "Loss" = "#B2182B",
      "No change" = "#E6C85C",
      "Gain" = "#4DAF4A"
    )
  ) +
  labs(
    x = NULL,
    y = "Change in suitable area relative to current (%)",
    colour = NULL,
    shape = "Scenario"
  ) +
  envalysis::theme_publish() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

ggsave(filename = 'figures/whole_win_los.png', 
       dpi=600,
       width = 9, height = 8, scale=1.5, units = 'cm')

evals <- read_csv("outputs/csv/plants_ensemble_eval.csv")

glimpse(evals)
glimpse(areas_models)

evals_clean <- evals %>% 
  mutate(
    species = str_remove(full.name, "_.*") %>% 
      str_replace_all("\\.", "_"),
    metric.eval = str_to_lower(metric.eval),
    filtered.by = str_to_lower(filtered.by)
  ) %>% 
  filter(filtered.by == "tss") %>% 
  select(
    species,
    metric.eval,
    cutoff,
    sensitivity,
    specificity,
    calibration
  ) %>% 
  mutate(species = str_replace(species, '_', ' ') %>% 
           str_to_sentence())

evals_clean %>% 
  filter(metric.eval == 'tss') %>% 
  distinct(species, sensitivity, 
           specificity, calibration) %>% 
  write_csv('outputs/summarise_metrics_eval_plants.csv')

write_csv(evals_clean, 'outputs/sup_mat_eval_metrics_plants.csv')

evals_wide <- evals_clean %>% 
  group_by(species, metric.eval) %>% 
  summarise(
    cutoff = mean(cutoff, na.rm = TRUE),
    sensitivity = mean(sensitivity, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    calibration = mean(calibration, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  pivot_wider(
    names_from = metric.eval,
    values_from = c(cutoff, sensitivity, specificity, calibration),
    names_glue = "{metric.eval}_{.value}"
  )

supp_long_table <- areas_models %>% 
  rename(
    p10_cutoff = cutoff
  ) %>% 
  left_join(evals_wide, by = "species") %>% 
  mutate(
    change_class_10pct = case_when(
      area_change_pct < -10 ~ "Loss",
      area_change_pct >  10 ~ "Gain",
      TRUE ~ "No change"
    )
  ) %>% 
  select(
    species,
    gcm,
    ssp,
    period,
    current_area_km2,
    future_area_km2,
    area_change_km2,
    area_change_pct,
    change_class_10pct,
    p10_cutoff,
    everything(),
    -path,
    -dir_path,
    -binary_file,
    -exists,
    -current_file,
    -occs_path,
    -kappa_cutoff,
    -roc_cutoff, -aucroc_cutoff, -aucroc_calibration,
    -kappa_sensitivity, -aucroc_specificity,
    -all_of(c('roc_sensitivity','aucroc_sensitivity', 'kappa_specificity', 
            'roc_specificity')) 
    
  ) %>% 
  arrange(species, ssp, period, gcm)

write_csv(supp_long_table, 'outputs/sup_area_change_plants.csv')

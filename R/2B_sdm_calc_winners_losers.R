library(tidyverse)
library(scales)
library(terra)
library(sf)
library(rnaturalearth)
library(scico)

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

areas_models <- read.csv("outputs/plants_projected_area_change_10p.csv")

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

scen_lab <- c(ssp245 = "Low-emissions\n(SSP245)",
              ssp370 = "Mid-emissions\n(SSP370)",
              ssp585 = "High-emissions\n(SSP585)")
per_lev  <- c("2015-2044", "2045-2074", "2075-2100")
gcms     <- c("CNRM-ESM2-1", "EC-Earth3-Veg", "INM-CM5-0", "MIROC6")

# species with no baseline area (NA / Inf change) cannot be plotted
areas_models %>% filter(!is.finite(area_change_pct)) %>% distinct(species)

# which combinations are repeated, how often, and for which species
dups <- areas_models %>% count(species, gcm, ssp, period) %>% filter(n > 1)
dups %>% count(copies = n)
dups %>% distinct(species) %>% pull(species)

# clean: exclude the endemic outlier, drop non-finite values, remove repeated rows
areas_clean <- areas_models %>%
  filter(species != "plagiobothrys_lithocaryus",
         is.finite(area_change_pct)) %>%
  distinct(species, gcm, ssp, period, .keep_all = TRUE)

# checks
n_distinct(areas_clean$species)                      
areas_clean %>% count(ssp, period) %>% distinct(n)   # expect 260 (65 x 4) in every cell
max(areas_clean$area_change_pct)                     # new maximum
set.seed(1)

pl <- areas_clean %>%
  mutate(
    ssp    = factor(ssp, levels = names(scen_lab)),
    period = factor(period, levels = per_lev),
    gcm    = factor(gcm, levels = gcms),
    # scenario bin + period offset + model offset + small jitter for the species
    x = as.integer(ssp) + (as.integer(period) - 2) * 0.30 +
      (as.integer(gcm) - 2.5) * 0.055 + runif(n(), -0.02, 0.02)
  )

# median across all species x models in each scenario x period cell
med <- pl %>%
  group_by(ssp, period) %>%
  summarise(med = median(area_change_pct), avg = mean(area_change_pct), .groups = "drop") %>%
  mutate(x = as.integer(ssp) + (as.integer(period) - 2) * 0.30)

p <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -10, ymax = 10,
           fill = "grey90", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_jitter(data = pl, colour = 'transparent',
             aes(x, area_change_pct, fill = period, shape = gcm),
             size = 2.2, alpha = 0.35) +
  geom_boxplot(data = pl, 
               aes(x, area_change_pct, fill = period, 
                   group=interaction(period, ssp)),
               outlier.size = 1.2,
               size = 0.3, alpha = 0.65) +
  geom_point(data = med,
             aes(x = x, y = avg),
             shape = 18, size=2.6, fill='black') +
  scale_x_continuous(breaks = 1:3, labels = unname(scen_lab),
                     expand = expansion(add = 0.5)) +
  scale_y_continuous(transform = transform_pseudo_log(sigma = 100),
                     breaks = c(-100, -50, -25, -10, 0, 10, 25, 50, 100, 250, 500)) +
  scale_fill_scico_d(palette = "batlow", begin = 0.1, end = 0.75, name = "Period",
                     labels = \(x) sub("-", "\u2013", x)) +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Global climate model") +
  guides(fill  = guide_legend(order = 1, override.aes = list(shape = 21, size = 3, alpha = 1)),
         shape = guide_legend(order = 2, override.aes = list(fill = "grey70", size = 3, alpha = 1))
  ) +
  labs(x = NULL, y = "Change in suitable area (%)") +
  theme_classic()+
  theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin())

p

ggsave("figures/fig3_plants_area_change_models.pdf", p,scale = 1.4,
       width = 16, height = 14, units = "cm", dpi = 300)
# species count behind the figure (expected: 300 / 4 = 75)
n_distinct(pl$species)

# what produces the extreme gains: baseline areas of the top gainers
areas_models %>%
  filter(is.finite(area_change_pct)) %>%
  arrange(desc(area_change_pct)) %>%
  select(species, gcm, ssp, period, current_area_km2, future_area_km2, area_change_pct) %>%
  head(10)

# projections above the axis maximum shown in the figure
ymax <- 600   # set to your axis limit
pl %>% filter(area_change_pct > ymax) %>% count(ssp, period)

# --- 1. Pivot to wide format ---------------------------------------------
# One row per species, columns = ssp x period combinations
wide_tbl <- area_change_plot %>%
  mutate(species = str_to_sentence(str_replace_all(species, "_", " "))) %>%
  pivot_wider(
    names_from  = c(ssp, period),
    values_from = area_change_pct,
    names_sep   = "_"
  ) %>%
  arrange(species)

# --- 2. Reorder columns explicitly (ssp245, ssp370, ssp585 x 3 periods each) ---
periods <- c("2015-2044", "2045-2074", "2075-2100")
ssps    <- c("ssp245", "ssp370", "ssp585")
ordered_cols <- as.vector(sapply(ssps, function(s) paste(s, periods, sep = "_")))

wide_tbl <- wide_tbl %>%
  select(species, all_of(ordered_cols))

# --- 3. Rename columns for display (just periods; SSP shown via header group) ---
display_names <- c("Species", rep(c("2015–2044", "2045–2074", "2075–2100"), times = 3))
colnames(wide_tbl) <- display_names

# --- 4. Build LaTeX table --------------------------------------------------
latex_code <- wide_tbl %>%
  kbl(
    format    = "latex",
    booktabs  = TRUE,
    longtable = TRUE,
    align     = c("l", rep("r", 9)),
    digits    = 1,
    caption   = "Projected percent change in habitat suitability area for each host-plant species across three SSP scenarios (SSP2-4.5, SSP3-7.0, SSP5-8.5) and three time periods.",
    label     = "tab:area_change_species"
  ) %>%
  kable_styling(latex_options = c("repeat_header", "striped"), font_size = 8) %>%
  add_header_above(c(" " = 1, "SSP2-4.5" = 3, "SSP3-7.0" = 3, "SSP5-8.5" = 3))

# --- 5. Write to file so you can copy-paste straight into your .tex -------
writeLines(latex_code, "plants_area_change_supp_table.tex")

# Print to console too
print(latex_code)

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
  filter(species != 'plagiobothrys_lithocaryus') %>% #removed outlier
  filter(ssp == 'ssp245', period == '2075-2100') %>% 
  select(species, gcm, period, area_change_pct) %>% 
  prep_counts('area_change_pct')

high_2100 <- areas_models %>%
  filter(species != 'plagiobothrys_lithocaryus') %>%  #removed outlier
  filter(ssp == 'ssp585', period == '2075-2100') %>% 
  select(species, gcm, period, area_change_pct) %>% 
  prep_counts('area_change_pct')

d <- areas_models %>%  
  filter(species != 'plagiobothrys_lithocaryus') %>%  #removed outlier
  select(species, ssp, gcm, period, area_change_pct) %>% 
  mutate(change_bin = cut(area_change_pct, breaks = bin_breaks,
                          labels = bin_labels, ordered_result = TRUE)) %>%
  mutate(change_bin = forcats::fct_expand(change_bin, bin_labels))

areas_models %>% count(species, ssp, gcm, period) %>% count(n)
areas_models %>%
  count(species, ssp, gcm, period, name = "n_rows") %>%
  filter(n_rows > 1) %>%
  distinct(species, n_rows)

d %>% filter(is.na(change_bin)) %>%
  summarise(n = n(), n_missing = sum(is.na(area_change_pct)),
            lo = min(area_change_pct, na.rm = TRUE),
            hi = max(area_change_pct, na.rm = TRUE))
range(bin_breaks)

sp <- d %>%
  distinct(species, ssp, gcm, period, .keep_all = TRUE) %>%   # after confirming duplicates
  filter(!is.na(area_change_pct)) %>%
  group_by(species, ssp, period) %>%
  summarise(med = median(area_change_pct),
            lo  = min(area_change_pct),
            hi  = max(area_change_pct), .groups = "drop")
thr <- 10   # % treated as no meaningful change

tally <- sp %>%
  mutate(dir = case_when(med <= -thr ~ "loss",
                         med >=  thr ~ "gain",
                         TRUE        ~ "stable")) %>%
  count(period, ssp, dir) %>%
  tidyr::pivot_wider(names_from = dir, values_from = n, values_fill = 0) %>%
  mutate(label = sprintf("%d loss · %d stable · %d gain", loss, stable, gain))
tally <- tally %>%
  left_join(sp %>% group_by(period, ssp) %>%
              summarise(disagree = sum(lo < 0 & hi > 0), .groups = "drop"),
            by = c("period", "ssp"))

scen_lab <- c(ssp245 = "Low-emissions\n(SSP245)",
              ssp370 = "Mid-emissions\n(SSP370)",
              ssp585 = "High-emissions\n(SSP585)")
lab <- tally %>% filter(period == '2075-2100') %>%
  mutate(txt = paste0(scen_lab, "\n", label)) %>%
  select(ssp, txt) %>% tibble::deframe()

sp %>% 
  filter(period == '2075-2100') %>% 
  ggplot(aes(x = med, y = ssp, colour = med)) +
  annotate('rect', xmin=-10, xmax=10, ymin=-Inf, ymax=Inf, 
           fill='gray80',  linewidth = NA, inherit.aes = F,alpha = 0.35) +
  geom_vline(xintercept = 0, linewidth = .3, colour = "grey40", linetype = 2) +
  geom_linerange(aes(xmin = lo, xmax = hi),
                 position = ggbeeswarm::position_quasirandom(orientation = "y", width = .3),
                 alpha = .25, linewidth = .5) +
  ggbeeswarm::geom_quasirandom(orientation = "y", size = 3, alpha = .7, width = .3) +
  # facet_wrap(~ period, ncol = 1) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 10),
                     breaks = c(-100, -50, -25, 0, 25, 50, 100, 200, 400)) +
  scale_colour_gradientn(colours = scico(n = 15, direction = 1, palette = "oleron"), 
                         rescaler = ~ scales::rescale_mid(.x, mid = 18)
                         )+
    
    # values=scico(n=6, begin =0, end=1, alpha = 1, categorical = T,
                                 # direction = -1, palette = 'tokyo')[c(6,1,5)])+
  scale_y_discrete(labels = lab) +
  labs(x = "Change in suitable area (%)", y = NULL) +
  envalysis::theme_publish(base_size = 14) +
  coord_flip()+
  theme(
        legend.position = 'none',      # top-right inside panel
        legend.justification = c(1, 1))

ggsave(filename = 'figures/win_los2.pdf', 
       dpi=600,
       width = 13, height = 10, scale=1.6, units = 'cm')

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

# --- 1. Overall percentages (across all species/ssp/period combos) --------
overall_pct <- area_change_plot_clean %>%
  count(change_class) %>%
  mutate(pct = 100 * n / sum(n))

overall_pct

# --- 2. Percentages grouped by period --------------------------------------
period_pct <- area_change_plot_clean %>%
  count(period, change_class) %>%
  group_by(period) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

period_pct |> group_by(change_class) |> mutate(mean_pct = mean(pct))

# --- 3. (Optional) wide format for easy reporting in text/table -----------
period_pct_wide <- period_pct %>%
  select(period, change_class, pct) %>%
  tidyr::pivot_wider(names_from = change_class, values_from = pct)

period_pct_wide

pal   <- scico(15, palette = "oleron")
brks  <- c(-100, -50, -25, 0, 25, 50, 100, 200, 400)
ylims <- range(c(sp$lo, sp$hi), na.rm = TRUE)   # shared across panels

make_panel <- function(p, first = FALSE) {
  
  lab <- tally %>%
    filter(period == p) %>%
    mutate(txt = paste0(ssp, "\n", label)) %>%
    select(ssp, txt) %>%
    tibble::deframe()
  
  ggplot(filter(sp, period == p), aes(x = ssp, y = med, colour = med)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -10, ymax = 10,
             fill = "grey80", colour = NA, alpha = .35) +
    geom_hline(yintercept = 0, linewidth = .3, colour = "grey40", linetype = 2) +
    geom_linerange(aes(ymin = lo, ymax = hi),
                   position = ggbeeswarm::position_quasirandom(orientation = "x", width = .3),
                   alpha = .25, linewidth = .4) +
    ggbeeswarm::geom_quasirandom(orientation = "x", width = .3, size = 2.8, alpha = .7) +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 10), breaks = brks) +
    scale_x_discrete(labels = lab) +
    scale_colour_gradientn(colours = pal,
                           limits = c(-100, 100), oob = squish,
                           rescaler = ~ rescale_mid(.x, mid = 0)) +
    coord_cartesian(ylim = ylims) +
    labs(x = NULL, title = p,
         y = if (first) "Change in suitable area (%)" else NULL) +
    envalysis::theme_publish(base_size = 12) +
    theme(legend.position = "none",
          plot.title   = element_text(face = "bold", size = rel(1)),
          axis.text.y  = if (first) element_text() else element_blank(),
          axis.ticks.y = if (first) element_line() else element_blank())
}

periods <- sort(unique(sp$period))
panels  <- purrr::imap(periods, ~ make_panel(.x, first = TRUE))
library(patchwork)

fig_s <- wrap_plots(panels, nrow = 3) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
fig_s

ggsave(filename = 'figures/whole_win_los.png', 
       dpi=600,
       width = 14, height = 18, scale=1.5, units = 'cm')

# ggplot(area_change_plot_clean, aes(x = period, y = area_change_pct)) +
#   geom_hline(
#     yintercept = 0,
#     linewidth = 0.4,
#     linetype = "dashed",
#     colour = "grey40"
#   ) +
#   geom_hline(
#     yintercept = c(-10, 10),
#     linewidth = 0.25,
#     linetype = "dotted",
#     colour = "grey60"
#   ) +
#   geom_violin(
#     aes(group = period),
#     width = 0.65,
#     fill = "grey90",
#     colour = "grey55",
#     linewidth = 0.3,
#     alpha = 0.45,
#     trim = FALSE
#   ) +
#   geom_jitter(
#     aes(shape = ssp, colour = change_class),
#     width = 0.12,
#     height = 0,
#     alpha = 0.65,
#     size = 1.9,
#     stroke = 0.6
#   ) +
#   stat_summary(
#     fun = median,
#     geom = "point",
#     shape = 95,
#     size = 8,
#     linewidth = 0.9,
#     colour = "black"
#   ) +
#   scale_y_continuous(
#     labels = label_percent(scale = 1)
#   ) +
#   scale_colour_manual(
#     values = c(
#       "Loss" = "#B2182B",
#       "No change" = "#E6C85C",
#       "Gain" = "#4DAF4A"
#     )
#   ) +
#   labs(
#     x = NULL,
#     y = "Change in suitable area relative to current (%)",
#     colour = NULL,
#     shape = "Scenario"
#   ) +
#   envalysis::theme_publish() +
#   theme(
#     legend.position = "bottom",
#     legend.box = "vertical"
#   )



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

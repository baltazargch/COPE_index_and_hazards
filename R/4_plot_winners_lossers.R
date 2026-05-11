library(tidyverse)
library(scales)
library(terra)
library(sf)

# library(gt)
evals <- read_csv('outputs/csv/plants_ensemble_eval.csv') %>% 
  filter(metric.eval == 'TSS') %>% 
  mutate(species = str_remove(full.name, '_.*') %>% 
           str_replace_all('\\.', '_'), .before=1) %>% 
  distinct(species, cutoff)

csv_mods <- read_csv('outputs/models/plants_projection_table.csv')

csv_mods <- csv_mods %>% 
  mutate(
    dir_path = str_remove(path, basename(path)),
    binary_file = path %>% basename() %>% 
      sub("^([^_]+)_(.+\\.tif)$", "\\1_TSSbinay_\\2", .) %>% 
      str_c(dir_path, .), 
    exists = file.exists(binary_file), 
    currrent = paste0('outputs/models/', species, '_NA_cur_ensemble.tif')
  ) %>% 
  left_join(evals, by='species')


# World basemap via rnaturalearth (land polygons)
wrld_org <- rnaturalearth::ne_states(country = 'united states of america', 
                                     returnclass = "sf")

ca <- wrld_org %>% filter(name_en == 'California')


curr <- rast(csv_mods$currrent[1]) %>% mean() %>% crop(ca, mask=T) 

plot( curr >= csv_mods$cutoff[1])


plot( curr >= csv_mods$cutoff[1])


res <-  list.files('outputs/models/plants_projections/', 
                   '.tif$', 
                  recursive = T, full.names = T)

prres <- read_csv('outputs/models/maxent/log_models_tunning.csv')


# res %>% 
#   select(current_area:last_col()) %>%
# mutate(
#   across(!matches('current|species'), \(x) (x - current_area)/current_area*100)
# ) %>% view
  

out <- res %>% 
  select(current_area:last_col()) %>% 
  mutate(
    across(!matches('current|species'), \(x) (x - current_area)/current_area*100)
  )

scenarios <- grep("^\\d{4}-\\d{4}_ssp\\d+$", names(res), value = TRUE)


# summary table: (current_area - scenario) for each scenario
summ_diff <- map_dfr(scenarios, function(scn) {
  diffs <- (res[[scn]] - res$current_area)/res$current_area*100
  sd <- sd(diffs, na.rm = T)
  
  s <- summary(diffs)  # Min., 1st Qu., Median, Mean, 3rd Qu., Max (and NA's)
  tibble(
    scenario = scn,
    Min      = unname(s["Min."]),
    Q1       = unname(s["1st Qu."]),
    Median   = unname(s["Median"]),
    Mean     = unname(s["Mean"]),
    sd = sd,
    Q3       = unname(s["3rd Qu."]),
    Max      = unname(s["Max."]),
    NAs      = if ("NA's" %in% names(s)) unname(s["NA's"]) else 0L
  )
})

summ_diff %>% select(-NAs)

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

ssp585_2015 <- prep_counts(out, "2075-2100_ssp245") 
ssp585_2075 <- prep_counts(out, "2075-2100_ssp585")


ssp585a <- ssp585_2015 %>% 
  filter(!is.na(Case)) %>% 
  ggplot(aes(x = change_bin, y = n_species, fill = Case)) +
  geom_col(position = "dodge",  alpha=0.8) +
  # geom_vline(xintercept = 5, linetype = "dashed", color = "gray12") +
  envalysis::theme_publish(base_size = 14) +
  scale_fill_manual(values=c('#B2182B', '#E6C85C', '#4DAF4A'),)+
  # ggsci::scale_fill_npg(alpha = 0.75)+
  labs(x = "% Change bin", y = "Number of species",  title = "Low-emissions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ssp585b <- ssp585_2075 %>% 
  filter(!is.na(Case)) %>% 
  ggplot(aes(x = change_bin, y = n_species, fill = Case)) +
  geom_col(position = "dodge",  alpha=0.8) +
  # geom_vline(xintercept = 5, linetype = "dashed", color = "gray12") +
  envalysis::theme_publish(base_size = 14) +
  scale_fill_manual(values=c('#B2182B', '#E6C85C', '#4DAF4A'),)+
  # ggsci::scale_fill_npg(alpha = 0.75)+
  labs(x = "% Change bin", y = "Number of species",  title = "High-emissions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

rbind(ssp585_2015 %>% mutate(Scenario = 'Low emissions', alpha='a') ,
      ssp585_2075 %>% mutate(Scenario = 'High emissions', alpha='b')) %>% 
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

terra::ext(terra::rast("outputs/monthly_climates/tasmin_2081_2100_ssp585.tif"))
res$species %>% unique() %>% paste(collapse = ', ')

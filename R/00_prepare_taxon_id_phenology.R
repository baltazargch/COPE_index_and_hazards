library(tidyverse)

plants <- read_csv('outputs/records/all_species_records_and_native.csv') %>% 
  rename_with(janitor::make_clean_names) %>% 
  filter(ca_native == 'yes')

plants_clean <- plants %>% 
  mutate(species = case_when(
    species == 'Mahonia aquifolium' ~ 'Berberis aquifolium', 
    .default = species))

cch2_taxonid <- read_csv('inputs/plants_phenology/Native Plants of California_1758044999.csv')

cch2_in <- cch2_taxonid %>% filter(ScientificName  %in% plants_clean$species)

plants_clean$species[!plants_clean$species  %in% cch2_taxonid$ScientificName] %>% 
  paste(., collapse = ', ')

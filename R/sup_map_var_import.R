library(tidyverse)
library(biomod2)

models <- list.files(
  "outputs/models/",
  pattern = "ensemble_model.rds",
  recursive = TRUE,
  full.names = TRUE
)

get_ensemble_varimp <- function(model_path) {
  
  bm_ens <- readRDS(model_path)
  
  species <- bm_ens@sp.name %>%
    str_replace_all("\\.", " ") %>%
    str_to_sentence()
  
  ens_varimp <- biomod2::bm_PlotVarImpBoxplot(
    bm.out   = bm_ens,
    group.by = c("expl.var", "algo", "full.name"),
    do.plot  = FALSE
  )
  
  ens_varimp$tab %>%
    mutate(
      species = species,
      model_path = model_path,
      .before = 1
    ) %>%
    rename(variable = expl.var) %>%
    group_by(species, variable) %>%
    summarise(
      mean_importance = mean(var.imp, na.rm = TRUE) * 100,
      sd_importance   = sd(var.imp, na.rm = TRUE) * 100,
      n_models        = n(),
      .groups = "drop"
    ) %>%
    mutate(
      mean_importance = round(mean_importance, 2),
      sd_importance   = round(sd_importance, 2)
    ) %>%
    arrange(species, desc(mean_importance))
}

varimp_long <- map_dfr(models, get_ensemble_varimp)

var_order <- varimp_long %>%
  mutate(has_data = mean_importance > 0) %>%
  group_by(variable) %>%
  summarise(n_species = sum(has_data)) %>%
  arrange(desc(n_species)) %>%
  pull(variable)

varimp_long_plot <- varimp_long %>%
  mutate(variable = factor(variable, levels = var_order))

ggplot(varimp_long_plot, aes(x = variable, y = species, fill = mean_importance)) +
  geom_tile() +
  geom_point(aes(size = sd_importance),shape=21, fill='white', color = "black") +
  scale_size(name = 'S.D.')+
  scale_fill_viridis_c(name = "Importance (%)") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    axis.text.y = element_text(face = 'italic')
  ) + 
  labs(x = 'Variable', y='Species')

ggsave(filename = 'suplemental_variable_importance.svg', 
       height = 20, width = 15, scale = 0.7)

varimp_complete <- varimp_long %>%
  complete(species, variable, fill = list(mean_importance = 0))

var_summary <- varimp_complete %>%
  group_by(variable) %>%
  summarise(
    overall_mean = mean(mean_importance, na.rm = TRUE),
    sd_between_species = sd(mean_importance, na.rm = TRUE),
    n_species = n_distinct(species),
    n_used = sum(mean_importance > 0),
    prop_used = n_used / n_species,
    .groups = "drop"
  ) %>%
  arrange(desc(overall_mean))


top_vars <- var_summary %>%
  slice_max(overall_mean, n = 5)

top_names <- top_vars$variable

glue::glue("
Across species, variable importance was dominated by a subset of climatic predictors. 
On average, the most influential variables were {top_names[1]}, {top_names[2]}, and {top_names[3]}, 
with mean importance values of {round(top_vars$overall_mean[1],1)}%, 
{round(top_vars$overall_mean[2],1)}%, and {round(top_vars$overall_mean[3],1)}%, respectively. 
These variables were consistently represented across species (present in 
{round(top_vars$prop_used[1]*100,0)}–{round(top_vars$prop_used[3]*100,0)}% of models), 
indicating shared climatic constraints. 

In contrast, other predictors showed lower mean importance and greater variability among species 
(sd = {round(mean(var_summary$sd_between_species, na.rm = TRUE),1)}%), suggesting more 
species-specific responses to environmental conditions.
")

res <- read_csv("outputs/final_maxent/csv/final_models_metrics.csv")
res$bin_path
res %>% select(!any_of('pred_path', 'model_rds', 'bin_path'))
colnames(res)

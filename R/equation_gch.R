library(dplyr)
library(tidyr)
library(ggplot2)



cope <- function(sf, rf, rc, k, v = 0.5, delta = 1, h = 3){
  
  # g_rc <- rf / (rc + h)
  g_rc <- rc / (rc + h)
  # g_rc <- rc + rf / (rc + rf + h)
  
  wr <- 1 + v * g_rc * tanh(k * log((rf + delta)/(rc + delta)))
  # wr <- 1 + v * tanh(k * log((rf + delta)/(rc + delta)))
  
  cope <- sf * wr
  return(cope)
}

# current richness values to compare
rc_vals <- seq(1, 300, by= 1)

# full continuum of future richness values
rf_vals <- 1:300

df <- expand.grid(
  rc = rc_vals,
  rf = rf_vals,
  k = c(0.2, 1), # 0.5, 0.7, 1, 2, 3, 5),
  sf = seq(0,1,by=0.01)) %>%
  as_tibble() %>%
  filter(rf >= 0) %>%
  mutate(
    case = paste0("Current richness = ", rc),
    delta_r = rf - rc,
    ratio = (rf + 1) / (rc + 1),
    wr = 1 + 0.5 * (rc / (rc + 5)) * tanh(k * log((rf + 1)/(rc + 1))),
    cope = cope(sf, rf, rc, k),
    delta_cope = cope - sf,
    cope_raw = cope_raw(sf, rf, rc, k),
    delta_cope_raw = cope_raw - sf
  )


df_main <- df %>%
  filter(sf == 0.6) 

label_df <- df_main %>%
  distinct(case, sf) %>%
  mutate(
    x = 5,
    y = 0.48,
    label = glue::glue("Sf = {sf}"))

library(tidyverse)


library(dplyr)
library(forcats)
library(ggplot2)
library(glue)

sel_var <- 'cope'

df_plot <- df_main %>%
  mutate(case = fct_reorder(case, rc)) %>%
  # filter(k %in% c(0.5, 2, 5))
  filter(k %in% c(0.2, 0.5, 0.7, 1))

df_labels <- df_plot %>%
  group_by(case, k, sf) %>%
  slice_max(order_by = delta_r, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = round( .data[[sel_var]], 2))

df_labels_min <- df_plot %>%
  group_by(case, k, sf) %>%
  slice_min(order_by = delta_r, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = round( .data[[sel_var]], 2))

# df_plot %>%
#   ggplot(aes(x = delta_r, y = !!sym(sel_var), colour = factor(k))) +
#   # geom_hline(yintercept = 1, linewidth = 0.2, linetype = "dashed") +
#   geom_hline(aes(yintercept = sf), linetype = 4, colour = "gray18") +
#   geom_vline(xintercept = 0, linewidth = 0.2, linetype = "dashed") +
#   geom_line(linewidth = 0.75) +
#   geom_label(data = df_labels,
#     aes(x = max(delta_r) - 2, y = sf, label = sf),
#     colour='gray8',
#     hjust = -0.2,
#     size = 3,
#     show.legend = FALSE
#   ) +
#   geom_label(
#     data = df_labels,
#     aes(label = label),
#     hjust = -0.2,
#     size = 3,
#     show.legend = FALSE
#   ) +
#   geom_label(
#     data = df_labels_min,
#     aes(label = label),
#     hjust = 1.5,
#     size = 3,
#     show.legend = FALSE
#   ) +
#   ggsci::scale_color_aaas() +
#   scale_x_continuous(expand = expansion(mult = c(0.25, 0.18)))+
#   facet_wrap(~case, scales = "free_x") +
#   labs(
#     x = "Change in resource richness (Rf - Rc)",
#     y = "COPE",
#     title = glue("Future suitability Sf:{unique(df_labels$sf)}\n", 
#                  "{sel_var}"),
#     colour = expression(kappa)
#   ) +
#   coord_cartesian(clip = "off") +
#   theme_light(base_size = 12) +
#   theme(
#     plot.margin = margin(5.5, 25, 5.5, 5.5)
#   ) 
df$k %>% range()
df %>% 
  filter(k %in% c(0.2, 0.8, 1), 
         # between(ratio, -0.1,0.1)
         ) %>% 
  group_by(k) %>% 
  sample_n(1e05) %>% 
  ggplot(aes(x = sf, y = delta_cope, 
             # size = log(rc),
             colour=log(ratio)))+
  geom_point(alpha=0.5, size=1.2) + 
  scale_colour_gradient2(low   = 'darkred', 
                       mid  = 'gray88', 
                       high ='forestgreen') + 
  facet_wrap(~k)


df_main$k %>% unique()
kin <- c(0.2, 0.5, 0.7, 1)
rcin <- sample(rc_vals, 10)
df_plot %>%
  filter(
    rc  %in% rcin,
    k %in% kin
  ) %>% 
  ggplot(aes(x = delta_r, y = !!sym(sel_var), group=case,
             colour = rc)) +
  geom_hline(aes(yintercept = sf), linetype = 4, colour = "gray18") +
  geom_vline(xintercept = 0, linewidth = 0.2, linetype = "dashed") +
  geom_line(linewidth = 0.75, na.rm = T) +
  geom_label(data = df_labels,
             aes(x = max(delta_r) - 100, y = sf, label = sf),
             colour='gray8',
             hjust = -0.2,
             size = 3,
             show.legend = FALSE
  ) +
  # geom_label(
  #   data = df_labels %>%  filter(
  #     rc  %in% rcin, 
  #     k %in% kin
  #   ) ,
  #   aes(label = label),
  #   hjust = -0.2,
  #   size = 3,
  #   show.legend = FALSE
  # ) +
  # geom_label(
  #   data = df_labels_min %>%  filter(
  #     rc  %in% rcin, 
  #     k %in% kin
  #   ) ,
  #   aes(label = label),
  #   hjust = 1.5,
  #   size = 3,
  #   show.legend = FALSE
  # ) +
  # ggsci::scale_color_aaas(na.translate = F) +
  facet_wrap(~k, scales = "free_x") +
  labs(
    x = "Change in resource richness (Rf - Rc)",
    y = "COPE",
    title = glue("Future suitability Sf:{unique(df_labels$sf)}\n", 
                 "{sel_var}"),
    colour = 'Case'#expression(kappa)
  ) +
  coord_cartesian(clip = "off") +
  theme_light(base_size = 12) +
  theme(
    plot.margin = margin(5.5, 25, 5.5, 5.5)
  ) 


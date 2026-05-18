library(tidyverse)
library(patchwork)
library(scales)

#-----------------------------
# Functions
#-----------------------------
wr_fun <- function(rf, rc, k, v = 0.5, delta = 1, h = 3){
  g_rc <- rc / (rc + h)
  1 + v * g_rc * tanh(k * log((rf + delta) / (rc + delta)))
}

cope_fun <- function(sf, rf, rc, k, v = 0.5, delta = 1, h = 3){
  sf * wr_fun(rf = rf, rc = rc, k = k, v = v, delta = delta, h = h)
}

#-----------------------------
# Palettes
#-----------------------------
pal_purple <- c(
  "#F3E8FF", "#E0AAFF", "#C77DFF", "#9D4EDD", "#7B2CBF", "#5A189A", "#240046"
)

pal_sf <- c(
  "Low suitability (Sf = 0.2)"  = "#E0AAFF",
  "Mid suitability (Sf = 0.5)"  = "#9D4EDD",
  "High suitability (Sf = 0.8)" = "#5A189A"
)

pal_rc <- c(
  "Low baseline richness (Rc = 10)"   = "#E0AAFF",
  "Moderate baseline richness (Rc = 50)"  = "#C77DFF",
  "High baseline richness (Rc = 150)" = "#7B2CBF",
  "Very high baseline richness (Rc = 300)" = "#240046"
)

# parameter values
k_vals <- c(0.1, 0.5, 0.9)   # if you want only endpoints, use c(0, 1)
h_vals <- 1:10
rc_cases <- c(10, 50, 150, 300)
ratio_vals <- seq(0.1, 3, length.out = 250)

wr_lines <- crossing(
  rc = rc_cases,
  ratio = ratio_vals,
  k = k_vals,
  h = h_vals
) %>%
  mutate(
    rf = ratio * (rc + 1) - 1,
    rf = pmax(rf, 0),
    delta_r = rf - rc,
    wr = wr_fun(rf = rf, rc = rc, k = k, h = h),
    # sf_lab = case_when(
    #   sf == 0.2 ~ "Low suitability (Sf = 0.2)",
    #   sf == 0.5 ~ "Mid suitability (Sf = 0.5)",
    #   sf == 0.8 ~ "High suitability (Sf = 0.8)"
    # ),
    # 
    rc_lab = case_when(
      rc == 10  ~ "Low baseline richness (Rc = 10)",
      rc == 50  ~ "Moderate baseline richness (Rc = 50)",
      rc == 150 ~ "High baseline richness (Rc = 150)",
      rc == 300 ~ "Very high baseline richness (Rc = 300)"
    ),
    rc_lab = factor(rc_lab, levels = names(pal_rc)),
    k_lab = factor(paste0("k = ", k), levels = paste0("k = ", k_vals))
  )

wr_lines_sum <- wr_lines %>%
  group_by(rc_lab, ratio, delta_r, k_lab) %>%
  summarise(
    wr_mean = mean(wr, na.rm = TRUE),
    wr_min = min(wr, na.rm = TRUE),
    wr_max = max(wr, na.rm = TRUE),
    .groups = "drop"
  )

p_wr_robust <- ggplot(wr_lines_sum, aes(x = ratio, y = wr_mean, colour = rc_lab, fill = rc_lab)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = wr_min, ymax = wr_max), alpha = 0.35, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_grid(k_lab~rc_lab, 
             labeller = labeller(
               rc_lab = label_wrap_gen(width = 20),
               k_lab  = label_wrap_gen(width = 15))) +
  scale_colour_manual(values = pal_rc) +
  scale_fill_manual(values = pal_rc) +
  # scale_x_continuous(
  #   breaks = c(0.25, 0.5, 1, 2, 3),
  #   labels = c("0.25", "0.5", "1", "2", "3")
  # ) +
  labs(
    x = expression((R[f] + delta) / (R[c] + delta)),
    y = expression(omega[r]),
    colour = NULL,
    fill = NULL,
    title = expression("Effect of richness change on " * omega[r]),
    subtitle = "Lines show the mean response; ribbons show the range across h = 1 to 10"
  ) +
  envalysis::theme_publish() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

p_wr_robust

h_surface <- c(1, 5, 10)
k_surface <- c(0.1, 0.5, 0.9)

wr_surface <- crossing(
  rc = seq(1, 300, by = 4),
  delta_r = seq(-250, 250, by = 4),
  k = k_surface,
  h = h_surface
) %>%
  mutate(
    rf = pmax(rc + delta_r, 0),
    wr = wr_fun(rf = rf, rc = rc, k = k, h = h),
    h_lab = factor(paste0("h = ", h), levels = paste0("h = ", h_surface)),
    k_lab = factor(paste0("k = ", k), levels = paste0("k = ", k_surface))
  )

p_wr_surface <- ggplot(wr_surface, aes(x = delta_r, y = rc, fill = wr)) +
  geom_raster() +
  facet_grid(k_lab ~ h_lab) +
  scale_fill_gradient2(
    low = "#F3E8FF",
    mid = "white",
    high = "#5A189A",
    midpoint = 1,
    name = expression(omega[r])
  ) +
  labs(
    x = expression(R[f] - R[c]),
    y = expression(R[c]),
    title = expression("Parameter surface of " * omega[r]),
    subtitle = "Values above 1 indicate a positive resource effect; values below 1 indicate a negative effect"
  ) +
  envalysis::theme_publish() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold")
  )

p_wr_surface

sf_vals <- c(0.2, 0.5, 0.8)
k_vals <- c(0.1, 0.5, 0.9)

cope_lines <- crossing(
  rc = rc_cases,
  ratio = ratio_vals,
  sf = sf_vals,
  k = k_vals,
  h = h_vals
) %>%
  mutate(
    rf = ratio * (rc + 1) - 1,
    rf = pmax(rf, 0),
    delta_r = rf - rc,
    cope = cope_fun(sf = sf, rf = rf, rc = rc, k = k, h = h),
    delta_cope = cope - sf,
    
    sf_lab = case_when(
      sf == 0.2 ~ "Low suitability (Sf = 0.2)",
      sf == 0.5 ~ "Mid suitability (Sf = 0.5)",
      sf == 0.8 ~ "High suitability (Sf = 0.8)"
    ),
    
    rc_lab = case_when(
      rc == 10  ~ "Low baseline richness (Rc = 10)",
      rc == 50  ~ "Moderate baseline richness (Rc = 50)",
      rc == 150 ~ "High baseline richness (Rc = 150)",
      rc == 300 ~ "Very high baseline richness (Rc = 300)"
    ),
    
    sf_lab = factor(sf_lab, levels = names(pal_sf)),
    rc_lab = factor(rc_lab, levels = names(pal_rc)),
    k_lab = factor(paste0("k = ", k), levels = paste0("k = ", k_vals))
  )


cope_lines_sum <- cope_lines %>%
  group_by(rc_lab, sf_lab, ratio, delta_r, k_lab) %>%
  summarise(
    delta_cope_mean = mean(delta_cope, na.rm = TRUE),
    delta_cope_min = min(delta_cope, na.rm = TRUE),
    delta_cope_max = max(delta_cope, na.rm = TRUE),
    .groups = "drop"
  )

p_cope <- ggplot(cope_lines_sum, aes(x = ratio, y = delta_cope_mean, 
                                     colour = sf_lab, fill = sf_lab)) +
  geom_ribbon(aes(ymin = delta_cope_min, ymax = delta_cope_max), 
              alpha = 0.15, colour = NA) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.9) +
  facet_grid(k_lab~ rc_lab, 
             labeller = labeller(
               rc_lab = label_wrap_gen(width = 20),
               k_lab  = label_wrap_gen(width = 15)
             )) +
  scale_colour_manual(values = pal_sf) +
  scale_fill_manual(values = pal_sf) +
  # scale_x_continuous(
  #   breaks = c(0.25, 0.5, 1, 2, 3),
  #   labels = c("0.25", "0.5", "1", "2", "3")
  # ) +
  labs(
    x = expression((R[f] + delta) / (R[c] + delta)),
    y = expression(COPE - S[f]),
    colour = NULL,
    fill = NULL,
    title = "Effect of resource change on COPE",
    subtitle = "Lines show the mean response; ribbons show the range across h = 1 to 10"
  ) +
  envalysis::theme_publish() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )


p_wr_robust
ggsave("figures/cope_parameter_k_richenss cases.png", p_wr_robust, 
       width = 12, height = 10, dpi = 300,scale = 0.8)

p_wr_surface
ggsave("figures/cope_parameter_wr_surface.png", p_wr_surface, 
       width = 12, height = 10, dpi = 300,scale = 0.8)

p_cope

ggsave("figures/cope_parameter_exploration_deltaCOPE.png", p_cope,
       width = 12, height = 10, dpi = 300, scale = 0.8)

# =========================================================
# Parameter space exploration for:
#   C = Sf * omega_r
#   omega_r = 1 + v * tanh(kappa * log((Rf + delta)/(Rc + delta)))
# =========================================================

library(tidyverse)
library(patchwork)
library(scales)

# ---------------------------------------------------------
# 1. Parameters
# ---------------------------------------------------------

# Baseline parameters
Rc    <- 100
delta <- 1
v     <- 0.5
h <- 3

# Focus on realistic changes in resource richness:
# Rf / Rc from ~0.1 to 2.0
ratio_seq <- seq(0.25, 1.75, length.out = 300)

# Translate ratio into Rf assuming Rc = 1
Rf_seq <- ratio_seq * Rc

# Explore kappa broadly, with emphasis on interpretable values
kappa_vals <- seq(0,1,by=0.1)

# Sf range
Sf_seq <- seq(0.1, 1.0, length.out = 200)

# Selected resource-ratio scenarios to show impact on C
ratio_selected <- c(0.2, 0.4, 0.8, 1.0, 1.25, 1.5, 1.75)

# ---------------------------------------------------------
# 2. Helper function
# ---------------------------------------------------------

omega_fun <- function(Rf, Rc = 1, delta = 0.05, v = 0.6, 
                      h = 3, kappa = 1) {
  1 + v * (Rc/(Rc + h)) * tanh(kappa * log((Rf + delta) / (Rc + delta)))
}

# ---------------------------------------------------------
# 3. Data for panel 1:
#    omega_r as a function of resource ratio, by kappa
# ---------------------------------------------------------

df_omega_curve <- crossing(
  ratio = ratio_seq,
  kappa = kappa_vals
) %>%
  mutate(
    Rf = ratio * Rc,
    omega_r = omega_fun(Rf = Rf, Rc = Rc, h = h, delta = delta, v = v, kappa = kappa),
    change_pct = (ratio - 1) * 100,
    kappa_f = factor(kappa, levels = kappa_vals)
  )

# ---------------------------------------------------------
# 4. Data for panel 2:
#    heatmap of omega_r across ratio x kappa
# ---------------------------------------------------------

df_omega_heat <- crossing(
  ratio = seq(0.25, 1.75, length.out = 240),
  kappa = seq(0, 1, by=0.1)
) %>%
  mutate(
    Rf = ratio * Rc,
    omega_r = omega_fun(Rf = Rf, Rc = Rc, delta = delta, v = v, kappa = kappa),
    change_pct = (ratio - 1) * 100
  )

# ---------------------------------------------------------
# 5. Data for panel 3:
#    effect on C across Sf for selected ratios and kappa
# ---------------------------------------------------------

df_C_curve <- crossing(
  Sf = Sf_seq,
  ratio = ratio_selected,
  kappa = kappa_vals
) %>%
  mutate(
    Rf = ratio * Rc,
    omega_r = omega_fun(Rf = Rf, Rc = Rc, delta = delta, v = v, kappa = kappa),
    C = Sf * omega_r,
    ratio_lab = paste0("Rf/Rc = ", ratio),
    kappa_f = factor(kappa, levels = kappa_vals)
  )

# ---------------------------------------------------------
# 6. Palette
# ---------------------------------------------------------
# Custom palette: blues + purples + fuchsia
#kappa_vals <- c(0.5, 1.5, 2.5, 4, 6)
pal_kappa <- c(
  "0"   = "#5DA5DA",  # soft blue
  "0.3"   = "#3B82F6",  # vivid blue
  "0.5"   = "#7C3AED",  # purple
  "0.7"     = "#C026D3",  # fuchsia-purple
  "1"     = "#EC4899"   # pink/fuchsia
)

# Continuous palette for heatmap
pal_heat <- c(
  "#EAF2FF",  # very light blue
  "#C7DBFF",
  "#9EC5FE",
  "#6FA8FF",
  "#5B8FF9",
  "#7C3AED",
  "#9D2EEA",
  "#C026D3",
  "#E043B6",
  "#F472B6"
)

# Ratio palette for selected scenarios
#ratio_selected <- c(0.2, 0.4, 0.8, 1.0, 1.25, 1.5, 2.0)

pal_ratio <- c(
  "Rf/Rc = 0.2"  = "#9EC5FE",  # strong blue
  "Rf/Rc = 0.4"  = "#5B8FF9",  # lighter blue
  "Rf/Rc = 0.8"  = "#7C3AED",  # purple
  "Rf/Rc = 1"    = "#9333EA",  # violet (center reference)
  "Rf/Rc = 1.25" = "#C026D3",  # fuchsia-purple
  "Rf/Rc = 1.5"  = "#E11D8A",  # strong fuchsia
  "Rf/Rc = 1.75" = "#F472B6"   # light magenta
)

# ---------------------------------------------------------
# 7. Plot theme
# ---------------------------------------------------------

theme_param <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# ---------------------------------------------------------
# 8. Panel 1: omega_r vs resource ratio
# ---------------------------------------------------------

p1 <- ggplot(df_omega_curve, aes(x = ratio, y = omega_r, color = kappa_f)) +
  geom_line(linewidth = 1.3) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey50", linewidth = 0.5) +
  scale_color_manual(values = pal_kappa, name = expression(kappa)) +
  scale_x_continuous(
    breaks = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75),
    labels = label_number(accuracy = 0.01)
  ) +
  scale_y_continuous(breaks = seq(0.5,1.5,by=0.1))+
  labs(
    title = expression("A. Response of " * omega[r] * " to resource richness ratio"),
    subtitle = bquote(v == .(v) * ", " ~ delta == .(delta) * ", " ~ Rc == .(Rc)),
    x = expression("Resource richness ratio (" * R[f] / R[c] * ")"),
    y = expression(omega[r])
  ) +
  theme_param
# p1
# ---------------------------------------------------------
# 9. Panel 2: heatmap of omega_r across ratio and kappa
# ---------------------------------------------------------
midpoint <- 1

p2 <- ggplot(df_omega_heat, aes(x = ratio, y = kappa, fill = omega_r)) +
  geom_raster(interpolate = TRUE) +
  
  geom_contour(
    aes(z = omega_r),
    color = "white",
    alpha = 0.45,
    bins = 8,
    linewidth = 0.35
  ) +
  
  scale_fill_gradient2(
    low = "#9EC5FE",      # blue
    mid = "#7C3AED",      # purple (neutral)
    high = "#EC4899",     # fuchsia
    # midpoint = midpoint,
    name = expression(omega[r])
  ) +
  
  scale_x_continuous(
    breaks = c(0.25,0.5,0.75,1,1.25,1.5,1.75),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  
  labs(
    title = expression("B. Parameter surface of " * omega[r]),
    subtitle = expression("Higher " * kappa * " sharpens the transition around " * R[f]/R[c] == 1),
    x = expression("Resource richness ratio (" * R[f] / R[c] * ")"),
    y = expression(kappa)
  ) +
  
  theme_param
p2
# ---------------------------------------------------------
# 10. Panel 3: effect on C across Sf
# ---------------------------------------------------------
df_C_curve$ratio %>% unique()

p3 <- ggplot(df_C_curve, aes(x = Sf, y = C, color = ratio_lab)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ kappa_f, nrow = 1, 
             labeller = labeller(kappa_f = function(x) paste("k = ", x))) +
  scale_color_manual(values = pal_ratio, name = expression(R[f] / R[c])) +
  labs(
    title = expression("C. Effect of " * omega[r] * " on multiplier " * C * " across " * S[f]),
    subtitle = "Same Sf can yield very different C depending on resource-ratio context and kappa",
    x = expression("Suitability " * S[f]),
    y = expression(C)
  ) +
  theme_param +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = NA)
  )

# ---------------------------------------------------------
# 11. Combine with patchwork
# ---------------------------------------------------------

final_plot <- (p1 | p2) / p3 +
  plot_annotation(
    title = "Parameter space exploration of the resource-response multiplier",
    subtitle = paste0(
      "Exploration of how kappa controls the sensitivity of omega_r to changes in resource richness, ",
      "and how this propagates to C through Sf"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 11, color = "grey30")
    )
  )

# Print
final_plot

# ---------------------------------------------------------
# 12. Save
# ---------------------------------------------------------

ggsave(
  filename = "figures/parameter_space_kappa_exploration.png",
  plot = final_plot,
  width = 13,
  height = 15,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "figures/parameter_space_kappa_exploration.svg",
  plot = final_plot,
  width = 13,
  height = 15,
  dpi = 600,
  bg = "white"
)

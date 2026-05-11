library(terra)
library(tidyverse)
set.seed(4525)

suit <- rast('outputs/models/CA_osli_current_ensemble.tif')

rt0 <-  rast('outputs/current_native_plant_richness_TSS.tif')
rt1 <-  sample(1:300, 35)

p_d <- 1
drp <- (rt1 +  p_d) / ( rt0 + p_d)

dr <- rt1 - rt0
cols <- ifelse(dr < 0, 'loss', 'gain')
cols[ between(dr, -5, 5) ] <- 'stable'


log_r <- log(drp)

plot(sort(log_r))

k = 2.9

for(k in seq(0.5,3,0.2) ){
  dp <-  tanh(round(k, 2)*log_r)
  
  plot(sort(dp), sort(dr))
  abline(h = 0, lty=2)
  title(paste0('k= ', round(k, 2)))
}


plot ( sort(atanh( (runif(100, 0,1) ))) )


h = 0.5
w = 1 + h * dp

plot( sort( dp ))


# jaccard -----------------------------------------------------------------
Ct0 <- lapply(1:35, \(c) sample(letters, 20))
Ct1 <- lapply(Ct0, \(c) sample(c( sample(setdiff(letters, c), 6), sample(c, 18)), 20))

# Function for Jaccard index
jaccard <- function(x, y) {
  inter <- length(intersect(x, y))
  union <- length(union(x, y))
  inter / union
}

sim <- runif(35, 0.3, 0.9)

plot(sort(sim))

max(sim); min(sim)

dp_j <- 2*sim-1
cbind(sim, dp_j)

db <- tibble(
  ini_rich = rt0, 
  end_rich = rt1,
  delta_prop = drp, 
  delta_abs = dr,
  log_dr = log_r, 
  suit = suit,
  # w_j = 1 + h * dp_j,
  dstar = dp, 
  w = w, 
  suit_w = (suit * w - min(suit * w, na.rm = TRUE)) / (max(suit * w, na.rm = TRUE) - min(suit * w, na.rm = TRUE)),
  col = cols
)

db %>% 
  ggplot(aes(suit, suit * w, size = abs(dr), colour = col))+
  geom_segment(aes(x = suit, y = suit, yend = suit * w), 
               inherit.aes = F, linewidth = 0.1)+
  geom_point() +
  geom_point(aes(suit, suit), inherit.aes = F)+
  ggsci::scale_color_aaas(alpha=0.8) + 
  theme_light()

db$w_j
db %>% 
  ggplot(aes(suit, suit * sim, size = abs(dr), colour = sim < 0.5 ))+
  geom_segment(aes(x = suit, y = suit, yend = suit* sim), 
               inherit.aes = F, linewidth = 0.1)+
  geom_point() +
  geom_point(aes(suit, suit), inherit.aes = F)+
  ggsci::scale_color_aaas(alpha=0.8) + 
  theme_light()


j <- seq(0, 1, 0.1)
dstar <- 2*j - 1
cbind(j, dstar)

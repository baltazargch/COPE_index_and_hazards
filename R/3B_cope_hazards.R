library(terra)
library(tidyverse)
library(networkD3)


# helper functions -------------------------------------------------------
cope_index <- function(sf, rf, rc, k, v = 0.5, delta = 1, h = 3){
  
  
  g_rc <- rc / (rc + h)
  
  wr <- 1 + v * g_rc * tanh(k * log((rf + delta)/(rc + delta)))
  
  cope <- sf * wr
  return(cope)
}

# helper: 4 quantiles (0–3)
cope_to_quants <- function(scop_r) {
  q <- quantile(scop_r[], na.rm = TRUE, probs = c(0, .25, .5, .75, 1))
  classify(
    scop_r,
    rcl = matrix(
      c(q[1], q[2], 0,
        q[2], q[3], 1,
        q[3], q[4], 2,
        q[4], q[5], 3),
      ncol = 3, byrow = TRUE
    ),
    right = TRUE, include.lowest = TRUE
  )
}

# ensure outside CA is NA (mask), inside CA 0/1
haz_bin_from_targets <- function(scope, ca, period) {
  # period = 'early'
  # scope = q
  
  fires <- rast(
    paste0('inputs/fires/', period, 'fires_rasters_intens_duration.tif')
  )
  
  fires.duration <- project(fires$duration, scope, method='near')
  fires.intensity <- project(fires$intenisty,scope, method='bilinear')
  
  v.fd <- values(fires.duration)
  v.fs <- values(fires.intensity) >= 2
  
  v.fd[ !v.fs ] <- NaN
  v.fd[ v.fs ] <- 1
  
  target.fires <- rast(fires.intensity, vals = v.fd) %>% mask(scope)
  
  # debugging
  # plot(target.fires)
  
  ### heatwaves ---------------------------------------------------------------
  hw <- list.files('inputs/heatwaves/', 'RCP85.*.rds', full.names = T) %>% 
    map(readRDS) %>% map(~.x[[ period ]]) %>% rast() %>% mean() %>% round(0)
  
  hw <- project(hw, scope, method='near')
  
  hw.masked <- crop(hw, ca, mask=T)
  
  terra::NAflag(hw.masked) <- 0
  
  v.hw <- values(hw.masked)
  vq <- quantile(v.hw, na.rm=T)[[4]]
  
  target.hw <- hw.masked >= vq
  # terra::NAflag(target.hw) <- 0
  
  target.hw <- target.hw %>% mask(scope)
  
  # hw_target and fire_target should be logical or 0/1 rasters already on same grid
  h <- (max(target.hw) == 1) | target.fires == 1
  vh <- values(h)
  vh[ is.na(vh) ] <- 0
  values(h) <- vh
  h <-  h %>% mask(scope)
  names(h) <- period
  h
}

# helper: build combined category raster
# hazard_r: 0 = no hazard, 1 = hazard (binary); NA outside CA
# quant_r: 0..3
# category code: hazard*4 + quant  => 0..7
build_category <- function(quant_r, hazard_r) {
  stopifnot(nlyr(quant_r) == 1, nlyr(hazard_r) == 1)
  cat_r <- hazard_r * 4 + quant_r
  names(cat_r) <- "cat"
  cat_r
}

# nice labels for plotting
cat_labels_binary <- function() {
  tibble(
    hazard_bin = rep(c("no.hazard", "hazard"), each = 4),
    scop_q = rep(c("no", "low", "mid", "high"), times = 2),
    cat = rep(0:7, 1)
  ) %>%
    mutate(label = paste(hazard_bin, scop_q, sep = ": "))
}

# cell.sizes <- cellSize(cope_ssps[[1]], unit="km")

transition_area <- function(cat_from, cat_to, area_r) {
  inter <- cat_from * 100 + cat_to
  names(inter) <- "inter"
  
  # sum area by inter-code
  z <- zonal(area_r, inter, fun = "sum", na.rm = TRUE)
  colnames(z) <- c("inter", "area_km2")
  
  z %>%
    as_tibble() %>%
    mutate(
      from = floor(inter / 100),
      to   = inter %% 100
    ) %>%
    select(from, to, area_km2)
}


make_sankey_df <- function(scop_list, area_r, periods = c("early","mid","late")) {
  # scop_list = cope_ssps
  # haz_bin_list = haz_bin
  # area_r = cell.sizes
  # periods = periods
  
  # 1) quantile class + combined category per period
  cats <- map2(periods,list(ca), \(p, c) {
    q <- cope_to_quants(scop_list[[p]])
    h <- haz_bin_from_targets(q, ca, p)
    vh <- values(h)
    
    vh[ is.na(vh) ] <- 0
    values(h) <- vh
    h <- crop(h, c, mask=T)
    build_category(q, h) 
  })
  
  names(cats) <- periods
  
  # 2) transitions for early->mid and mid->late
  flows_1 <- transition_area(cats[[periods[1]]], cats[[periods[2]]], area_r) %>%
    mutate(p_from = periods[1], p_to = periods[2])
  
  flows_2 <- transition_area(cats[[periods[2]]], cats[[periods[3]]], area_r) %>%
    mutate(p_from = periods[2], p_to = periods[3])
  
  flows <- bind_rows(flows_1, flows_2)
  
  # 3) create node names that include the period (so they’re separate columns in Sankey)
  lab <- cat_labels_binary()
  
  flows <- flows %>%
    left_join(lab, by = c("from" = "cat")) %>% rename(from_lab = label) %>%
    left_join(lab, by = c("to" = "cat"))   %>% rename(to_lab   = label) %>%
    mutate(
      source_name = paste0(p_from, " | ", from_lab),
      target_name = paste0(p_to,   " | ", to_lab)
    )
  
  nodes <- tibble(name = unique(c(flows$source_name, flows$target_name))) %>%
    mutate(id = row_number() - 1)
  
  links <- flows %>%
    left_join(nodes, by = c("source_name" = "name")) %>% rename(source = id) %>%
    left_join(nodes, by = c("target_name" = "name")) %>% rename(target = id) %>%
    select(source, target, value = area_km2)
  
  list(nodes = nodes, links = links)
}


# Load data and run -------------------------------------------------------

## calc scop_list ----------------------------------------------------------
ca <- rnaturaleplotca <- rnaturalearth::ne_states('United States of America')
ca <- ca %>% filter(name == 'California')

#this is current richness. Immutable
rc <- rast('outputs/baseline_floral_richness.tif') 

#files for future suitability
sf <- rast(list.files("outputs/models/projections/means/", 'ssp585_', full.names = T )) %>% 
  crop(ca, mask=T)/1000

#Calculate future richness
csv_mods <- read_csv('outputs/database_floral_models_10p_th.csv')

periods <- csv_mods$period %>% unique()

names(sf) <- periods

rf_all <- map(periods, \(p){
  # p=periods[1] #debug
  
  # future richness patterns for period
  ssp585_p_plants <- csv_mods %>% 
    filter(ssp=='ssp585', period == p) %>% 
    distinct(species, path) %>% 
    left_join(
      plants_bin_10th %>% map_dfr(\(x) tibble(species = x$species, cutoff = x$cutoff)), by='species'
    )
  
  plan(multisession, workers=18)
  p_rast <- future_map(unique(ssp585_p_plants$species), \(x){
    # x=ssp585_plants$species[1]
    sp_db <- ssp585_p_plants %>% filter(species == x) 
    
    high_emissions_map <- (rast(sp_db$path) %>% mean()/1000) %>% crop(ca, mask=T)
    
    ca_high_emissions_bin <- high_emissions_map >= sp_db$cutoff[1]
    wrap(ca_high_emissions_bin)
  })
  
  plan(sequential)
  
  period_richn <- map(p_rast, unwrap) %>% rast() %>% sum(na.rm = T) 
  names(period_richn) <- p
  period_richn
})

rf_all <- rast(rf_all)


cope_ssps <-  c(
  early = cope_index(sf[[1]], rf_all[[1]], rc, k = 0.5, v = 0.5, delta = 1, h = 3),
  mid = cope_index(sf[[2]], rf_all[[2]], rc, k = 0.5, v = 0.5, delta = 1, h = 3),
  late = cope_index(sf[[3]], rf_all[[3]], rc, k = 0.5, v = 0.5, delta = 1, h = 3)
)
       

names(cope_ssps) <- c('early', 'mid', 'late')

# debugging
# map(cope_ssps, plot)

cell.sizes <- cellSize(cope_ssps[[1]], unit='km')
## hazards -----------------------------------------------------------------
periods <- c('early', 'mid', 'late')

cats <- map2(periods,list(ca), \(p, c) {
  # p <- periods[1]
  # c = ca
  # 
  q <- cope_to_quants(cope_ssps[[p]])
  vh <- values(q)
  
  vh[ vh == 0 ] <- NA
  values(q) <- vh

  h <- haz_bin_from_targets(q, c, p)
  
  cat <- build_category(q, h) 
  
})

lab <- cat_labels_binary()

map(1:3, \(i){
  areas.haz <- zonal(cell.sizes, cats[[i]], sum, unit='km') %>% 
    left_join(lab, by='cat') 
  
  areas.haz %>% 
    group_by(hazard_bin, scop_q) %>% 
    summarise(area = sum(area),
              prop = area/sum(areas.haz$area) * 100)
  
})

### fires -------------------------------------------------------------------
# sanity check
# plot(target.fires)
# plot(target.hw)


sank <- make_sankey_df(cope_ssps, cell.sizes, periods)


nodes_filt <- sank$nodes %>%
  filter(!str_detect(name, ": no$")) %>%   # drop "no"
  mutate(
    new_id = row_number() - 1,
    period = str_extract(name, "^(early|mid|late)")  # node group
  )

links_filt <- sank$links %>%
  left_join(nodes_filt %>% select(id, new_id), by = c("source" = "id")) %>%
  rename(source_new = new_id) %>%
  left_join(nodes_filt %>% select(id, new_id), by = c("target" = "id")) %>%
  rename(target_new = new_id) %>%
  filter(!is.na(source_new), !is.na(target_new)) %>%
  transmute(
    source = as.integer(source_new),
    target = as.integer(target_new),
    value  = value
  )
colourScale <- '
d3.scaleOrdinal()
  .domain(["early","mid","late"])
  .range(["#E0E0E0","#A8A8A8","#5A5A5A"])
'
library(htmlwidgets)

# to make sanky in web sankey page
map_chr(1:nrow(links_filt), \(x) {
  paste0(
    nodes_filt$name[ nodes_filt$new_id == links_filt$source[x] ], 
    ' COPE [', round(links_filt$value[x], 2), '] ', 
    nodes_filt$name[ nodes_filt$new_id == links_filt$target[x] ], 
    ' COPE'
  ) 
}) %>% paste0(collapse = '\n') %>% cat


p <- networkD3::sankeyNetwork(
  Links  = links_filt,
  Nodes  = nodes_filt,
  Source = "source",
  Target = "target",
  Value  = "value",
  NodeID = "name",
  NodeGroup = "period",
  colourScale = colourScale,
  fontSize = 12,
  nodeWidth = 30
)
p

library(webshot2)

saveWidget(p, "outputs/sankey.html", selfcontained = TRUE)

# Try SVG (may work depending on setup)
webshot("outputs/sankey.html", file = "sankey.pdf", selector = "svg")

# Report results ----------------------------------------------------------
node_meta <- nodes_filt %>%
  separate(
    name,
    into = c("period", "hazard_conserv"),
    sep = " \\| ",
    remove = FALSE
  ) %>%
  separate(
    hazard_conserv,
    into = c("hazard", "conservation"),
    sep = ": ",
    remove = FALSE
  )

flows <- links_filt %>%
  left_join(node_meta, by = c("source" = "new_id")) %>%
  rename(
    period_from = period,
    hazard_from = hazard,
    conserv_from = conservation
  ) %>%
  left_join(node_meta, by = c("target" = "new_id")) %>%
  rename(
    period_to = period,
    hazard_to = hazard,
    conserv_to = conservation
  )

hazrds.id <- node_meta %>% filter(hazard == 'hazard') %>% pull(new_id)
# 
# links_filt %>% 
#   filter(hazard_from == 'hazard') %>% 
#   summarise(area_km2 = sum(value), .groups = "drop")

net_change <- flows %>%
  group_by(period_from, period_to, conserv_to) %>%
  summarise(area_km2 = sum(value), .groups = "drop")


transition_type <- flows %>%
  mutate(
    direction = case_when(
      conserv_to == conserv_from ~ "stable",
      conserv_to %in% c("mid", "high") & conserv_from == "low" ~ "improving",
      conserv_to == "high" & conserv_from == "mid" ~ "improving",
      conserv_to == "low" & conserv_from %in% c("mid", "high") ~ "degrading",
      conserv_to == "mid" & conserv_from == "high" ~ "degrading",
      TRUE ~ "other"
    )
  ) %>%
  group_by(period_to, direction) %>%
  summarise(area_km2 = sum(value), .groups = "drop")

hazard_effect <- flows %>%
  group_by(period_to, hazard_to, conserv_to) %>%
  summarise(area_km2 = sum(value), .groups = "drop")


haz_bin_list <- set_names(
  map(periods, \(.x) {
    h <- haz_bin_from_targets(scope = cope_ssps[[.x]], ca = ca, period = .x)
    vh <- values(h)
    
    vh[ is.na(vh) ] <- 0
    values(h) <- vh
    h <- crop(h, ca, mask=T)
  }),
  periods
)

# area per cell (km^2), using any template that matches your SCOP grid
area_km2_r <- cellSize(cope_ssps[[1]], unit = "km")

haz_areas <- imap_dfr(haz_bin_list, \(h, p) {
  # make sure grids align
  h2 <- project(h, area_km2_r, method = "near")
  
  # sum area grouped by hazard class (0/1)
  z <- zonal(area_km2_r, h2, fun = "sum", na.rm = TRUE)
  colnames(z) <- c("haz_bin", "area_km2")
  
  as_tibble(z) %>%
    mutate(
      period = p,
      exposure = if_else(haz_bin == 1, "hazard", "no.hazard")
    ) %>%
    select(period, exposure, area_km2)
})

haz_areas

haz_areas_pct <- haz_areas %>%
  group_by(period) %>%
  mutate(
    total_km2 = sum(area_km2),
    pct = 100 * area_km2 / total_km2
  ) %>%
  ungroup()

haz_areas_pct %>% 
  filter(period == 'late')

periods <- c("early", "mid", "late")

make_combined_map <- function(period) {
  scop_q <- cope_to_quants(cope_ssps[[period]]) |> crop(ca, mask = TRUE)
  scop_q[scop_q == 0] <- NA
  
  haz <- haz_bin_from_targets(scope = cope_ssps[[period]], ca = ca, period = period)
  haz <- project(haz, scop_q, method = "near")
  
  vh <- values(haz)
  
  vh[ is.na(vh) ] <- 0
  values(haz) <- vh
  haz <- crop(haz, ca, mask=T)
  
  
  comb <- haz * 10 + scop_q
  names(comb) <- period
  comb
}

class_levels <- c(
  "1"  = "No hazard: Low",
  "2"  = "No hazard: Mid",
  "3"  = "No hazard: High",
  "11" = "Hazard: Low",
  "12" = "Hazard: Mid",
  "13" = "Hazard: High"
)

comb_stack <- rast(map(periods, make_combined_map))
names(comb_stack) <- periods

vc <- values(comb_stack)
vc[is.na(vc)] <- 0

values(comb_stack) <- vc
comb_stack <- crop(comb_stack, ca, mask=T)

plot(comb_stack)

writeRaster(comb_stack %>% 
              disagg(fact=8, method='near') %>%
              mask(ca), 'outputs/cope_x_hazards.tif', overwrite=TRUE)

imap(cope_ssps,\(x,n) writeRaster(x%>% 
                                    disagg(fact=8, method='near') %>%
                                    mask(ca), paste0('outputs/', n, '_cope.tif'), overwrite=TRUE))

df_map <- as.data.frame(comb_stack, xy = TRUE, na.rm = TRUE) %>%
  pivot_longer(cols = all_of(periods), names_to = "period", values_to = "class") %>%
  mutate(
    class = factor(as.character(class), levels = names(class_levels), labels = unname(class_levels)),
    period = factor(period, levels = periods)
  )

scop_hazard_palette <- c(
  # no hazard (bright)
  "No hazard: Low"  = "#E5D8BD",  # light gray
  "No hazard: Mid"  = "#A6BD6D",  # mustard yellow
  "No hazard: High" = "#4DAF4A",  # medium green
  
  # hazard (dark)
  "Hazard: Low"     = "#E6C85C",  # dark gray
  "Hazard: Mid"     = "#F7C7C7",  # dark mustard
  "Hazard: High"    = "#B2182B"   # dark green
)
library(sf)
ggplot(df_map) +
  geom_raster(aes(x = x, y = y, fill = class)) +
  geom_sf(data = st_as_sf(ca), fill = NA, linewidth = 0.3) +
  # coord_equal() +
  facet_wrap(~ period, nrow = 1) +
  scale_fill_manual(
    values = scop_hazard_palette,
    name = "Conservation opportunity\nand hazard exposure"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


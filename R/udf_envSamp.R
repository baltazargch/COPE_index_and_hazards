library(terra)
library(tidyverse)

envSample <- function(coord, rasters, res, do.plot = TRUE) {
  
  # coord <- st_occs
  # rasters <- spPreds[[varsCor$selected.vars]]
  # res <- 1
  # 1. Extract environmental values at coordinates
  # coord is data.frame(lon, lat) or matrix
  
  # Extract environmental values
  vals <- terra::extract(rasters, coord) %>% as_tibble() %>% dplyr::select(-ID)
  
  n <- ncol(vals)
  
  # Precompute bin indices
  bins_df <- map2_dfc(1:n, res, function(i, r) {
    col <- vals[[i]]
    min_val <- min(col, na.rm = TRUE) - 1
    bin <- floor((col - min_val) / r)
    tibble(!!paste0("bin", i) := bin)
  })
  
  # Build groupID as combined bin indices
  bins_df <- bins_df %>%
    group_by(across(everything())) %>%
    mutate(groupID = cur_group_id()) %>%
    ungroup()
  
  # Add groupID to real_p
  real_p <- bind_cols(coord, vals, bins_df)
  
  # Keep 1 point per bin
  selected_points <- real_p %>%
    group_by(groupID) %>%
    slice(1) %>%  # take first point in bin
    ungroup()
  
  coord_filter <- selected_points %>% dplyr::select(lon=X, lat=Y)
  rows_mask <- coord$X %in% coord_filter$lon & coord$Y %in% coord_filter$lat
  sum(rows_mask)
  # Plotting
  if (do.plot) {
    par(mfrow = c(1, 2), mar = c(4, 4, 0, 0.5))
    
    plot(vals[[1]], vals[[2]], pch = 19, col = "grey50", 
         xlab = "Filter 1", ylab = "Filter 2")
    points(vals[rows_mask, ], pch = 20, col = "forestgreen", cex=1.1)
    
    plot(coord$X, coord$Y, pch = 19, col = "grey50")
    points(coord_filter$lon, coord_filter$lat, pch = 20, col = "forestgreen", cex=1.1)
    
    par(mfrow = c(1, 1))
  }
  
  # Return rows mask
  return(rows_mask)
}
library(terra)

terra_biovars <- function(prec, tmin, tmax) {
  
  stopifnot(
    inherits(prec, "SpatRaster"),
    inherits(tmin, "SpatRaster"),
    inherits(tmax, "SpatRaster"),
    nlyr(prec) == 12,
    nlyr(tmin) == 12,
    nlyr(tmax) == 12
  )
  
  # Select the value from a raster stack based on a layer index raster
  extract_by_index <- function(r_stack, index_rast) {
    mat <- values(r_stack)
    idx <- values(index_rast)[, 1]
    out <- rep(NA_real_, nrow(mat))
    good <- !is.na(idx)
    out[good] <- mat[cbind(which(good), idx[good])]
    setValues(rast(r_stack, nlyr = 1), out)
  }
  
  # Get mean temp of warmest/coldest quarter
  quarter_temp <- function(tavg, which = "max") {
    indices <- embed(c(1:12, 1:2), 3)[, 3:1]
    qmeans <- lapply(1:12, function(i) mean(tavg[[indices[i, ]]]))
    qstack <- rast(qmeans)
    extract_by_index(qstack, if (which == "max") which.max(qstack) else which.min(qstack))
  }
  
  # Get precipitation of wettest/driest quarter
  quarter_prec <- function(prec, which = "max") {
    indices <- embed(c(1:12, 1:2), 3)[, 3:1]
    qsums <- lapply(1:12, function(i) sum(prec[[indices[i, ]]]))
    qstack <- rast(qsums)
    extract_by_index(qstack, if (which == "max") which.max(qstack) else which.min(qstack))
  }
  
  # Get precipitation of warmest/coldest quarter (based on temperature)
  quarter_prec_by_temp <- function(prec, tavg, which = "max") {
    indices <- embed(c(1:12, 1:2), 3)[, 3:1]
    qtemp <- lapply(1:12, function(i) mean(tavg[[indices[i, ]]]))
    qtemp_stack <- rast(qtemp)
    
    qprec <- lapply(1:12, function(i) sum(prec[[indices[i, ]]]))
    qprec_stack <- rast(qprec)
    
    idx <- if (which == "max") which.max(qtemp_stack) else which.min(qtemp_stack)
    extract_by_index(qprec_stack, idx)
  }
  
  tavg <- (tmin + tmax) / 2
  
  # BIO1: Annual Mean Temperature
  bio1 <- mean(tavg)
  
  # BIO2: Mean Diurnal Range
  bio2 <- mean(tmax - tmin)
  
  # BIO5 and BIO6 for later use
  bio5 <- max(tmax)  # Max temperature of warmest month
  bio6 <- min(tmin)  # Min temperature of coldest month
  
  # BIO7: Temperature Annual Range
  bio7 <- bio5 - bio6
  
  # BIO3: Isothermality
  bio3_vals <- values(bio2) / values(bio7) * 100
  bio3_vals[!is.finite(bio3_vals)] <- NA
  bio3 <- setValues(rast(bio2), bio3_vals)
  
  # BIO4: Temperature Seasonality (SD * 100)
  bio4 <- app(tavg, fun = sd) * 100
  
  # BIO8: Mean Temp of Wettest Month
  bio8 <- extract_by_index(tavg, which.max(prec))
  
  # BIO9: Mean Temp of Driest Month
  bio9 <- extract_by_index(tavg, which.min(prec))
  
  # BIO10: Mean Temp of Warmest Quarter
  bio10 <- quarter_temp(tavg, "max")
  
  # BIO11: Mean Temp of Coldest Quarter
  bio11 <- quarter_temp(tavg, "min")
  
  # BIO12: Annual Precipitation
  bio12 <- sum(prec)
  
  # BIO13: Precipitation of Wettest Month
  bio13 <- max(prec)
  
  # BIO14: Precipitation of Driest Month
  bio14 <- min(prec)
  
  # BIO15: Precipitation Seasonality (CV)
  bio15 <- (app(prec, sd) / mean(prec)) * 100
  
  # BIO16: Precip of Wettest Quarter
  bio16 <- quarter_prec(prec, "max")
  
  # BIO17: Precip of Driest Quarter
  bio17 <- quarter_prec(prec, "min")
  
  # BIO18: Precip of Warmest Quarter
  bio18 <- quarter_prec_by_temp(prec, tavg, "max")
  
  # BIO19: Precip of Coldest Quarter
  bio19 <- quarter_prec_by_temp(prec, tavg, "min")
  
  # Stack all results
  rast(list(
    bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10,
    bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19
  )) |> `names<-`(paste0("bio", 1:19))
  
  # biovar_names <- c(
  #   bio1  = "ann_mean_temp",
  #   bio2  = "mean_diurnal_range",
  #   bio3  = "isothermality",
  #   bio4  = "temp_seasonality",
  #   bio5  = "max_temp_warmest_month",
  #   bio6  = "min_temp_coldest_month",
  #   bio7  = "temp_annual_range",
  #   bio8  = "mean_temp_wettest_month",
  #   bio9  = "mean_temp_driest_month",
  #   bio10 = "mean_temp_warmest_qtr",
  #   bio11 = "mean_temp_coldest_qtr",
  #   bio12 = "ann_precip",
  #   bio13 = "precip_wettest_month",
  #   bio14 = "precip_driest_month",
  #   bio15 = "precip_seasonality",
  #   bio16 = "precip_wettest_qtr",
  #   bio17 = "precip_driest_qtr",
  #   bio18 = "precip_warmest_qtr",
  #   bio19 = "precip_coldest_qtr"
  # )
}

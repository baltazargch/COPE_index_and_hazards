library(curl)
library(tidyverse)


dwlnd <- rbind(
  expand_grid(
    models = 'CNRM-ESM2-1', 
    ssp = 'ssp370', 
    veg = 'BAU',
    base = 'https://data.pyregence.org/wg4/CEC-Preliminary/wrf_bc/CNRM-ESM2-1_ssp370_r1i1p1f2/LANDIS-II/Veg-BAU/Replicate_101/', 
    case = paste0('CA_fire_severity_', 2001:2100, '.tif')
  ), 
  expand_grid(
    models = 'CNRM-ESM2-1', 
    ssp = 'ssp370', 
    veg = 'High Ambition',
    base = 'https://data.pyregence.org/wg4/CEC-Preliminary/wrf_bc/CNRM-ESM2-1_ssp370_r1i1p1f2/LANDIS-II/Veg-HighAmbition/Replicate_101/',
    case = paste0('CA_fire_severity_', 2001:2100, '.tif')
  ), 
  expand_grid(
    models = 'CNRM-ESM2-1', 
    ssp = 'ssp585', 
    veg = 'High Ambition',
    base = 'http://ungoliant.ucmerced.edu/data/CEC-Preliminary/loca2/CNRM-ESM2-1_ssp585_r1i1p1f2/LANDIS-II/Veg-HighAmbition/Replicate_211/',
    case = paste0('CA_fire_severity_', 2001:2100, '.tif')
  ), 
  expand_grid(
    models = 'EC-Earth3-Veg', 
    ssp = 'ssp370', 
    veg = 'BAU',
    base = 'https://data.pyregence.org/wg4/CEC-Preliminary/wrf_bc/EC-Earth3-Veg_ssp370_r1i1p1f1/LANDIS-II/Veg-BAU/Replicate_101/',
    case = paste0('CA_fire_severity_', 2001:2100, '.tif')
  ), 
  expand_grid(
    models = 'EC-Earth3-Veg', 
    ssp = 'ssp370', 
    veg = 'High Ambition',
    base = 'https://data.pyregence.org/wg4/CEC-Preliminary/wrf_bc/EC-Earth3-Veg_ssp370_r1i1p1f1/LANDIS-II/Veg-HighAmbition/Replicate_101/',
    case = paste0('CA_fire_severity_', 2001:2100, '.tif')
  ),
  expand_grid(
    models = 'INM-CM5-0', 
    ssp = 'ssp245', 
    veg = 'High Ambition',
    base = 'http://ungoliant.ucmerced.edu/data/CEC-Preliminary/loca2/INM-CM5-0_ssp245_r1i1p1f1/LANDIS-II/Veg-HighAmbition/Replicate_211/',
    case = paste0('CA_fire_severity_', 2001:2100, '.tif')
  )
  
)

walk(1:nrow(dwlnd), 
     \(x){
       dir.create(outdir <- paste0('inputs/fires/', 
                         dwlnd$models[x], '/',
                         dwlnd$ssp[x], '/',
                         dwlnd$veg[x]),
                         showWarnings = F, 
                         recursive = T)
       outfile <- paste0(outdir,'/', dwlnd$case[x])
       if(file.exists(outfile)) return(NULL)
       curl_download(str_c(dwlnd$base[x], dwlnd$case[x]),
                     destfile = outfile)
     })


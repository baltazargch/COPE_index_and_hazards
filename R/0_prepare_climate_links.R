library(tidyverse)

base <- 'https://cirrus.ucsd.edu/~pierce/LOCA2/NAmer/'

models <- c("INM-CM5-0", "EC-Earth3-Veg", "MIROC6", "CNRM-ESM2-1")

periods <- c('1950-2014', 
             '2015-2044', 
             '2045-2074', 
             '2075-2100')

ssp <- c('historical', 'ssp245', 'ssp370', 'ssp585')

vars <- c('pr', 'tasmax', 'tasmin')


suf1 <- '/0p0625deg/r1i1p1f1/'
suf1b <- '/0p0625deg/r1i1p1f2/'
suf2 <- '.LOCA_16thdeg_v20240915.monthly.nc'
suf3 <- '.LOCA_16thdeg_v20220413.monthly.nc'

suf4 <- 'r1i1p1f1'
suf4b <- 'r1i1p1f2'
g <- expand_grid(
  base, models, periods, ssp, vars, suf1, suf2, suf4
)
g_filtered <- g %>% 
  filter((ssp == "historical" & periods == "1950-2014") |
           (ssp != "historical" & periods != "1950-2014")) %>% 
  mutate(
    suf2 = ifelse(vars != 'pr', suf3, suf2), 
    suf1 = ifelse(models == 'CNRM-ESM2-1', suf1b, suf1), 
    suf4 = ifelse(models == 'CNRM-ESM2-1', suf4b, suf4), 
  ) %>% 
  mutate(
    path = glue('{base}{models}{suf1}{ssp}/{vars}/{vars}.{models}.{ssp}.{suf4}.{periods}{suf2}')
         )
  
g_filtered %>% pull(path) %>% 
  write_lines(., 'outputs/monthly_climates/NAclim/links_download.txt')

g_filtered %>% count(models)

library(callr)

writeLines('', 'rerun.txt')

rerun <- file.exists('rerun.txt')
while(rerun){
  rscript('R/2_SDM_maxent_plants.R')
  rerun <- file.exists('rerun.txt')
}

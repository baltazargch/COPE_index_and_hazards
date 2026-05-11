# ---- Helper: project one GCM/period/ssp and write EMmean layer --------------
project_one <- function(projOL, bm_ens, mods_sgl, mods_ens,n.cores, 
                        path, gcm, ssp, period, names_rast, sel_vars, 
                        out_dir) {
 
  env <- rast(path)
  names(env) <- names_rast
  env <- env[[sel_vars]]                    # subset to selected vars
  env <- rotate(env)
  env <- project(env, "epsg:4326")
  
  ens <- BIOMOD_EnsembleForecasting(
    bm.em         = bm_ens,
    proj.name           = paste(ssp, period, gcm, sep = "_"),
    new.env             = env,
    models.chosen = mods_ens,
    metric.binary = "TSS",
    nb.cpu        = n.cores
  )
  
  r <- unwrap(ens@proj.out@val)
  bin <- ens@proj.out@link
  
  
  # Prefer the EMmean layer; fall back to the first if not found
  idx <- grep("EMmean", names(r), fixed = TRUE)
  r_em <- if (length(idx)) r[[idx[1]]] else r[[1]]
  
  out_path <- file.path(out_dir, glue("{bm_ens@sp.name}_{gcm}_{ssp}_{period}_EMmean.tif"))
  writeRaster(r_em, out_path, overwrite = TRUE)
  
  ibin <- grep('TSSbin', bin)
  rbin <- rast(bin[ibin])
  writeRaster(rbin, 
              file.path(out_dir, glue("{bm_ens@sp.name}_TSSbinay_{gcm}_{ssp}_{period}_EMmean.tif")), 
              overwrite = TRUE)
  
  tibble(gcm = gcm, ssp = ssp, period = period, path = out_path)
}
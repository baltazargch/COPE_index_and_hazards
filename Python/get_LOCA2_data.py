# from Pierce, D. W., D. R. Cayan, D. R. Feldman, and M. D. Risser, 2023: 
# Future Increases in North American Extreme Precipitation in CMIP6 Downscaled with 
# LOCA. J. Hydrometeor., 24, 951–975, https://doi.org/10.1175/JHM-D-22-0194.1.
import os
import pandas as pd
from climakitae.core.data_interface import get_data
import xarray as xr

# Output folder
climate_dir = "inputs/climates"
os.makedirs(climate_dir, exist_ok=True)

# Target models and parameters
modelsIN = ["INM-CM5-0", "EC-Earth3-Veg", "MIROC6", "CNRM-ESM2-1"]

# Variables and scenarios
variables = {
    "tasmax": "Maximum air temperature at 2m",
    "tasmin": "Minimum air temperature at 2m",
    "pr": "Precipitation (total)"
}
scenarios = {
    "historical": "Historical Climate",
    "ssp245": "SSP 2-4.5",
    "ssp370": "SSP 3-7.0",
    "ssp585": "SSP 5-8.5"
}
# Define  periods
periods = [(2000, 2020), (2021, 2040), (2041, 2060), (2061, 2080), (2081, 2100)]

for scen_key, scen_name in scenarios.items():
    for start_year, end_year in periods:
        for var_short, var_label in variables.items():
            print(f"Processing: {var_short}, {scen_key}, {start_year}-{end_year}")
            outname = f"{var_short}_{start_year}_{end_year}_{scen_key}.nc"
            outpath = os.path.join(climate_dir, outname)
            
            if os.path.exists(outpath):
              print(f"⏭ Skipping existing file: {outpath}")
              continue
            
            
            units = "degC" if "tas" in var_short else "mm"

            data = get_data(
                variable=var_label,
                downscaling_method="Statistical",
                resolution="3 km",
                timescale="monthly",
                cached_area="CA",
                units=units,
                time_slice=(start_year, end_year),
                scenario=[scen_name]
            )

            # Filter simulations by model
            selected_sims = [
                sim for sim in data.coords["simulation"].values
                if any(m in sim for m in modelsIN)
            ]
            data = data.sel(simulation=selected_sims)

            # Save
           
            data.to_netcdf(outpath)
            print(f"✔ Saved: {outpath}")

# 
# 
# 
# # Get matching simulation names
# sim_names = pr.coords['simulation'].values
# selected_sims = [sim for sim in sim_names if any(model in sim for model in modelsIN)]
# 
# tasmax = tasmax.sel(simulation=selected_sims)
# tasmin = tasmin.sel(simulation=selected_sims)
# pr = pr.sel(simulation=selected_sims)
# 
# # Convert to datasets for indexing
# tasmax = tasmax.to_dataset(name="tasmax")
# tasmin = tasmin.to_dataset(name="tasmin")
# pr = pr.to_dataset(name="pr")
# 
# # Stack time into (year, month)
# time_index = pd.to_datetime(tasmax.time.values)
# 
# tasmax['year'] = ('time', time_index.year) 
# 
# tasmax['month'] = ('time', time_index.month)
# tasmin['year'] = ('time', time_index.year)
# tasmin['month'] = ('time', time_index.month)
# pr['year'] = ('time', time_index.year)
# pr['month'] = ('time', time_index.month)
# 
# tasmax = add_time_coords(tasmax)
# tasmin = add_time_coords(tasmin)
# pr = add_time_coords(pr)
# 
# # Group and compute per-year mean BIOCLIMs, then average across years
# bios = []
# for sim in tasmax.simulation.values:
#     print(f"Processing {sim}...")
#     sim_bio = []
#     for y in np.unique(tasmax.year):
#         # Extract 12-month slices
#         tx = tasmax.tasmax.sel(simulation=sim).where(tasmax.year == y).groupby('month').mean('time').values
#         tn = tasmin.tasmin.sel(simulation=sim).where(tasmin.year == y).groupby('month').mean('time').values
#         tm = (tx + tn) / 2.
#         pp = pr.pr.sel(simulation=sim).where(pr.year == y).groupby('month').mean('time').values
# 
#         # Check data
#         if tx.shape[1] != 12 or tn.shape[1] != 12 or tm.shape[1] != 12 or pp.shape[1] != 12:
#             print(f"Skipping {sim} {y} — incomplete year")
#             continue
#         
#   
#         
#         bio = biovars(pre=pp, tmp=tm, tmn=tn, tmx=tx)
#         sim_bio.append(bio)
#         
#         np.mean(pp.reshape(4, 3), axis=1)
# 
#     # Convert yearly biovars to numpy array and average over time
#     sim_bio = np.array(sim_bio)
#     sim_bio_mean = sim_bio.mean(axis=0)
# 
#     bios.append((sim, sim_bio_mean))
# 
# # Build Dataset
# bio_ids = [f"bio{i+1:02d}" for i in range(19)]
# bio_ds = xr.Dataset()
# 
# for i, bio_name in enumerate(bio_ids):
#     bio_ds[bio_name] = ("simulation", np.array([b[1][i] for b in bios]))
# 
# bio_ds["simulation"] = ("simulation", [b[0] for b in bios])


# gtsm3-era5-nrt
Code for GTSMv3.0 near-real time simulations with ERA5 on Snellius (for Sea Level Monitor)

# Step 1: Download data
1a ERA5
1b tides

# Step 2: Preprocess ERA5 to FM format
- fix lat coordinates
- change standard names
- merge daily files to yearly/monthly
-
# Step 3: Start GTSM model simulation
- use singularity containers
- template files

# Step 4: Postprocessing of raw data
- compute residual water levels
- fix attributes
- compute annual means
- plot data

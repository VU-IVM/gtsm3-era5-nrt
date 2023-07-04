import xarray as xr
import matplotlib.pyplot as plt

file_nc = "/gpfs/work1/0/einf3499/model_runs_extended/slr_tide_surge_runs/model_input_ERA5_1951/output/gtsm_fine_0000_his.nc"

data_xr = xr.open_dataset(file_nc)

max_wl = data_xr.waterlevel.max(dim='time')
min_wl = data_xr.waterlevel.min(dim='time')

fig,ax = plt.subplots()
pc = ax.scatter(x=max_wl.station_x_coordinate,y=max_wl.station_y_coordinate, s=0.5, c=max_wl, vmax=2)
fig.colorbar(pc, ax=ax)
fig.savefig('maxwl.png')



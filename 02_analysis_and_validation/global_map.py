import matplotlib as mpl
import cartopy as crt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.geoaxes import GeoAxes
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
def global_map(ax):
    # Projections
    crg = crt.crs.PlateCarree() # the one we have defined the data
    crgp = crt.crs.Robinson() # the one to plot the data
    # Identify european 
    #lonmin = -35; lonmax = 50; latmin = 20; latmax = 80
    # possibly use region mask
    # Plotting
    ax.set_global()
    ax.set_extent([-180, 180, -60, 70], crg)            
    ax.coastlines(resolution='10m', color='gray', linewidth=0.5, alpha=0.8,zorder=6)
    ax.add_feature(crt.feature.LAND.with_scale('10m'),facecolor='gray',zorder=4,alpha=0.20)
    ax.add_feature(crt.feature.LAKES.with_scale('10m'), facecolor='gray',zorder=5,alpha=0.05)        
    ax.add_feature(crt.feature.OCEAN.with_scale('10m'), edgecolor='face', facecolor='white')
    # Plot variable
    #bs=ax.scatter(x=ds['station_x_coordinate'],y=ds['station_y_coordinate'],
    #                     s=100,c=ds,transform=crg)
    # Format lat lon grid    
    gl = ax.gridlines(crs=crg, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.yline = gl.xlines = True
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.top_labels = gl.right_labels = False
    gl.ylabel_style = {'size': 12, 'color': 'gray'}
    gl.xlabel_style = {'rotation': 0, 'size': 12, 'color': 'gray'}
    gl.top_labels = gl.right_labels = False
    gl.ylocator = mpl.ticker.FixedLocator(np.arange(-90.,91.,20))
    gl.xlocator = mpl.ticker.FixedLocator(np.arange(-180.,181.,40))
    return ax
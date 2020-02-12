# EMS-python

A Python package to accurately plot the geographic data from a EMS-produced NetCDF file.


## Requirements

This package is based on the following dependencies :
- [matplotlib](https://matplotlib.org/)
- [netcdf4-python](http://unidata.github.io/netcdf4-python/netCDF4/index.html)
- [xarray](http://xarray.pydata.org/en/stable/), which allows to interact with the NetCDF files on every other level
- [cartopy](https://scitools.org.uk/cartopy/docs/latest/index.html)


## Installation

To install this package, install the required dependencies with any method you see fit then execute 
```
git clone https://github.com/MartinGraive/ems-python.git
cd ems-python/
python setup.py install
```


## Usage

Use the function `plot_map` to display a heatmap of the region with relevant data as described in the arguments.

Example:
```
import xarray as xr
import matplotlib.pyplot as plt
from emsplot import plot_map

# eReefs GBR4 hydrodynamic v2.0 model data
# produced by the CSIRO Coastal Environment modelling team
ds = xr.open_dataset('http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_v2/gbr4_simple_2019-01.nc')
# variable: temperature ('temp')
# on January 1, 2019 (time index 0)
# at elevation 50cm (k index 45)
plot_map(ds, 'temp', {'time': 0, 'k': 45}, cmap='inferno')

plt.show()
```

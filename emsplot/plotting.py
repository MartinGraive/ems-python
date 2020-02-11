import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

import warnings


def _latlon_autodetect(data, varname):
    var_dims = data[varname].dims
    lat, lon = [], []
    for k, i in data.variables.items():
        if i.attrs.get('coordinate_type') == 'latitude':
            if set(i.dims).issubset(var_dims):
                lat.append(k)
        if i.attrs.get('coordinate_type') == 'longitude':
            if set(i.dims).issubset(var_dims):
                lon.append(k)
    if len(lat) == 0 or len(lon) == 0:
        raise KeyError("Autodetection failed.")
    elif len(lat) > 1 or len(lon) > 1:
        print("Please choose adequate coordinates:")
        print(*["{}: {}".format(k, i) for k, i in enumerate(lat)], sep='\n')
        c1 = int(input(">> "))
        print(*["{}: {}".format(k, i) for k, i in enumerate(lon)], sep='\n')
        c2 = int(input(">> "))
        lat, lon = lat[c1], lon[c2]
    else:
        lat, lon = lat[0], lon[0]
    return lat, lon


def _title_creation(data, indexes):
    title = ""
    relevant_info = []
    for k, i in indexes.items():
        if isinstance(i, int):  # .isel()
            relevant_info.append(k)

    # search in .coords if it is defined
    for k, i in (data.coords.items() if len(data.coords) > 0 else data.items()):
        if len(i.dims) == 1 and i.dims[0] in relevant_info:  # fixme: hack
            kk = i.dims[0]
            val = i[indexes[kk]].data
            try:
                val = np.datetime_as_string(val, 'auto')
            except TypeError:  # workaround because of buggy lib
                pass
            unit = i.units if hasattr(i, 'units') else ''
            title += "{} = {}{}, ".format(k, val, unit)

    return title.rstrip(', ')


# def _plot_from_var(var, ax, proj, add_cbar, ji, **kwargs):
#     title = ""
#     latlon = []
#     if len(var.coords.variables) != 0:
#         for k, i in var.coords.items():
#             if i.size == 1:
#                 title += "{} = {}{}, ".format(k,
#                                             i.data if k != 'time' else np.datetime_as_string(i.data, 'auto'),
#                                             i.units if hasattr(i, 'units') else '')
#             else:
#                 latlon.append(k)
#     else:  # if the coordinates are not stored in "coordinates"
#         latlon = [k for k in var.dims]
#     if len(latlon) != 2:
#         raise IndexError('Wrong coordinates')
#     if ji:
#         latlon.reverse()
#     try:
#         im = ax.pcolor(var[latlon[0]], var[latlon[1]], var.data, transform=proj,
#                        **kwargs)
#     except TypeError:
#         im = ax.pcolor(var[latlon[1]], var[latlon[0]], var.data, transform=proj,
#                        **kwargs)
#
#     ax.set_title(title.rstrip(', '))
#
#     if add_cbar:
#         cbar = plt.colorbar(im, ax=ax)
#         cbar.set_label("{} ({})".format(var.long_name, var.attrs.get('units')))
#     return ax


# TODO: add a fill_value option for ill-defined arrays
# TODO: allow for .sel(**indexes)
# TODO: possibility to force projection
def plot_map(data, varname, indexes, ax=None, latlon=None,
         add_coastline=True, add_cbar=True, add_title=True, **kwargs):
    """Plot the map of the selected variable.

    Parameters
    ----------
    data : xarray.DataSet
        Dataset containing all variables
    varname : str
        Name of the variable to plot
    indexes : dict
        Dictionary of coordinates to set
        To take the whole range of a coordinate into account, do not mention it
    ax : cartopy.mpl.geoaxes.GeoAxesSubplot
        Object to draw on
        None to create it
    latlon : tuple or None
        Name given to the latitude and longitude variables
        None to autodetect
    add_coastline : bool
        Whether to add coastal outlines
    add_cbar : bool
        Whether to add a colorbar
    add_title : bool
        Whether to add a title
    kwargs
        matplotlib.axes.Axes.pcolor kwargs

    Returns
    -------
    cartopy.mpl.geoaxes.GeoAxesSubplot
        Object drawn on
    """
    if isinstance(latlon, (tuple, list)):
        lat, lon = latlon
    else:
        lat, lon = _latlon_autodetect(data, varname)
    lat, lon = data[lat], data[lon]

    proj_type = lat.attrs.get('projection')  # geographic_crs_name ?
    if proj_type == 'geographic':
        proj = ccrs.PlateCarree()
    else:
        warnings.warn("Could not check projection type, "
                      "defaults to rectangular", SyntaxWarning)
        proj = ccrs.PlateCarree()
    if ax is None:
        ax = plt.axes(projection=proj)

    try:
        var = data[varname].isel(**indexes)
    except ValueError as e:
        e.args += ("valid dimensions are {}".format(data[varname].dims),)
        raise

    posindexes = dict((k, indexes[k]) for k in lat.dims if k in indexes)
    im = ax.pcolor(lon.isel(**posindexes), lat.isel(**posindexes), var.data,
                   transform=proj, **kwargs)

    if add_cbar:
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("{} ({})".format(var.long_name, var.attrs.get('units')))

    if add_title:
        ax.set_title(_title_creation(data, indexes))

    if add_coastline:
        ax.coastlines(resolution="50m")  # enough for gbr4
    return ax


if __name__ == '__main__':
    """Tests
    """
    import xarray as xr
    # ds = xr.open_dataset('/home/martin/Bureau/AIMS/recom_test/outputs/out_std.nc')
    ds = xr.open_dataset('http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_v2/gbr4_simple_2019-01.nc')

    # plot_map(ds, 'temp', {'record': 6, 'k_centre': 21})
    plot_map(ds, 'temp', {'time': 0, 'k': 45}, cmap='inferno')

    plt.show()
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colorbar import make_axes, ColorbarBase

import warnings


# TODO: compatibility with netCDF4.Dataset (attributes have different names)

def _latlon_autodetect(data, varname):
    def condition(i, string):
        c1 = i.attrs.get('coordinate_type') == string
        if c1:
            return c1
        ln = i.attrs.get('long_name')
        if ln is None:
            return c1
        return string in ln.lower()

    var_dims = data[varname].dims
    lat, lon = [], []
    for k, i in data.variables.items():
        if condition(i, 'latitude'):
            if set(i.dims).issubset(var_dims):
                lat.append(k)
        if condition(i, 'longitude'):
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
        if isinstance(i, (int, np.int64)):  # .isel()
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


def _dim_coo_split(data, indexes):
    # if hasattr(data, 'coordinates'):
    #     coo_list = data.coordinates
    # else:
    #     coo_list = list(data.coords)
    dcoo, ddim = dict(), dict()
    for k, i in indexes.items():
        if k[0] == '+':
            dcoo[k[1:]] = i
        elif k in list(data.dims):
            ddim[k] = i
        else:
            dcoo[k] = i

    dcoo2 = dict()
    for k, i in dcoo.items():
        dk = data[k]
        # fixme: if more than one dim
        dcoo2[dk.dims[0]] = np.minimum(np.searchsorted(dk, i), len(dk) - 1)

    return {**dcoo2, **ddim}


# TODO: add a fill_value option for ill-defined arrays
# TODO: allow for .sel(indexes)
# TODO: interp and from true data (e.g. metres rather than index)
# TODO: possibility to force projection
# TODO: automatically choose the best resolution for coastlines
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
    try:
        if isinstance(latlon, (tuple, list)):
            lat, lon = latlon
        else:
            lat, lon = _latlon_autodetect(data, varname)
        lat, lon = data[lat], data[lon]

        indexes = _dim_coo_split(data, indexes)
        var = data[varname].isel(indexes)
    except KeyError as e:
        e.args += ("valid variable names are "
                   "{}".format(tuple(data.data_vars) + tuple(data.coords)),)
        raise
    except ValueError as e:
        e.args += ("valid dimensions are {}".format(data[varname].dims),)
        raise

    proj_type = lat.attrs.get('projection')  # geographic_crs_name ?
    if proj_type == 'geographic':
        proj = ccrs.PlateCarree()
    else:
        warnings.warn("Could not check projection type, "
                      "defaults to rectangular", SyntaxWarning)
        proj = ccrs.PlateCarree()
    if ax is None:
        ax = plt.axes(projection=proj)

    posindexes = (dict((k, indexes[k]) for k in lon.dims if k in indexes),
                  dict((k, indexes[k]) for k in lat.dims if k in indexes))
    lon, lat = lon.isel(posindexes[0]), lat.isel(posindexes[1])
    if len(lon.dims) == 1:  # we expect that lon and lat have the same n of dims
        lon = np.tile(lon, (lat.shape[0], 1))
        lat = np.tile(lat.T, (lon.shape[1], 1)).T
    im = ax.pcolor(lon, lat, var.data,
                   transform=proj, **kwargs)

    if add_cbar:
        cbar = plt.colorbar(im, ax=ax)
        try:
            name = var.long_name
        except AttributeError:
            name = var.name
        cbar.set_label("{} ({})".format(name, var.attrs.get('units')))

    if add_title:
        ax.set_title(_title_creation(data, indexes))

    if add_coastline:
        ax.coastlines(resolution="50m")  # enough for gbr4

    return ax


def animate_map(data, varname, indexes, anim_dim, anim_range=None, interval=300, ax=None,
                latlon=None, add_coastline=True, add_cbar=True, add_title=True, **kwargs):
    """Animate the map of the selected variable.

    Parameters
    ----------
    data : xarray.DataSet
        Dataset containing all variables
    varname : str
        Name of the variable to plot
    indexes : dict
        Dictionary of coordinates to set
        To take the whole range of a coordinate into account, do not mention it
    anim_dim : str
        Name of the coordinate to update
    anim_range : iterable, int, generator or None
        Range of the updated coordinate
        None to take whole range
    interval : int
        Delay between frames in milliseconds
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
    matplotlib.animation.FuncAnimation
        Animation object
    """
    try:
        if isinstance(latlon, (tuple, list)):
            lat, lon = latlon
        else:
            lat, lon = _latlon_autodetect(data, varname)
        lat, lon = data[lat], data[lon]

        indexes = _dim_coo_split(data, indexes)
        var = data[varname].isel(indexes)
    except KeyError as e:
        e.args += ("valid variable names are "
                   "{}".format(tuple(data.data_vars) + tuple(data.coords)),)
        raise
    except ValueError as e:
        e.args += ("valid dimensions are {}".format(data[varname].dims),)
        raise

    proj_type = lat.attrs.get('projection')  # geographic_crs_name ?
    if proj_type == 'geographic':
        proj = ccrs.PlateCarree()
    else:
        warnings.warn("Could not check projection type, "
                      "defaults to rectangular", SyntaxWarning)
        proj = ccrs.PlateCarree()
    if ax is None:
        ax = plt.axes(projection=proj)
    if add_cbar:
        cax, _ = make_axes(ax)

    posindexes = (dict((k, indexes[k]) for k in lon.dims if k in indexes),
                  dict((k, indexes[k]) for k in lat.dims if k in indexes))
    lon, lat = lon.isel(posindexes[0]), lat.isel(posindexes[1])
    if len(lon.dims) == 1:  # we expect that lon and lat have the same n of dims
        lon = np.tile(lon, (lat.shape[0], 1))
        lat = np.tile(lat.T, (lon.shape[1], 1)).T

    try:
        name = var.long_name
    except AttributeError:
        name = var.name

    def update(frame):
        arr = var.isel({anim_dim: frame}).data
        ax.pcolor(lon, lat, arr, transform=proj, **kwargs)
        if add_title:
            ax.set_title(_title_creation(data, {**indexes, anim_dim: frame}))
        if add_cbar:
            cax.clear()
            cbar = ColorbarBase(cax, boundaries=np.unique(arr[~np.isnan(arr)]).flatten(), **kwargs)
            cbar.set_label("{} ({})".format(name, var.attrs.get('units')))

    if isinstance(anim_range, slice):
        anim_range = np.arange(var.sizes[anim_dim])[anim_range]
    ani = FuncAnimation(ax.get_figure(), update, anim_range, interval=interval)

    if add_coastline:
        ax.coastlines(resolution="50m")  # enough for gbr4

    return ani


def _adapt_type(value, length):
    if isinstance(value, int):
        return [value]
    elif isinstance(value, slice):
        return np.arange(length)[value]
    else:
        return value


# TODO: time labels factorization
# TODO: diagonals
# TODO: interp and from true data (e.g. degrees rather than index)
def plot_ts(data, varname, indexes, ax=None, latlon=None, **kwargs):
    """Plot the time series of the selected variable at different locations.

    Parameters
    ----------
    data : xarray.DataSet
        Dataset containing all variables
    varname : str
        Name of the variable to plot
    indexes : dict
        Dictionary of coordinates to set
        Same sizes or use a single int to take all range
    ax : matplotlib.axes._subplots.AxesSubplot
        Object to draw on
        None to create it
    latlon : tuple or None
        Name given to the latitude and longitude variables
        None to autodetect
    kwargs
        matplotlib.axes.Axes.plot kwargs

    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot
        Object drawn on
    """
    try:
        if isinstance(latlon, (tuple, list)):
            lat, lon = latlon
        else:
            lat, lon = _latlon_autodetect(data, varname)
        lat, lon = data[lat], data[lon]
        posindexes = dict((k, indexes[k]) for k in lat.dims if k in indexes)

        indexes = _dim_coo_split(data, indexes)
        var = data[varname].isel(dict((k, indexes[k])
                                      for k in indexes if k not in posindexes))
    except KeyError as e:
        e.args += ("valid variable names are "
                   "{}".format(tuple(data.data_vars) + tuple(data.coords)),)
        raise
    except ValueError as e:
        e.args += ("valid dimensions are {}".format(data[varname].dims),)
        raise

    if ax is None:
        ax = plt.axes()

    plt.ylabel("{} ({})".format(var.long_name, var.attrs.get('units')))

    time_dim = [dim for dim in var.dims if dim not in lat.dims][0]
    for k, i in (data.coords.items() if len(data.coords) > 0 else data.items()):
        if len(i.dims) == 1 and i.dims[0] == time_dim:  # fixme: hack
            j = i[indexes[time_dim]]
            plt.xticks(np.arange(len(j.data)), np.datetime_as_string(j, 'auto'),
                       rotation=90)

    unitlat = lat.units if hasattr(lat, 'units') else ''
    unitlon = lon.units if hasattr(lon, 'units') else ''
    value_list = [_adapt_type(v, var.sizes[k]) for k, v in posindexes.items()]
    l0, l1 = len(value_list[0]), len(value_list[1])
    if l0 != l1:
        value_list = [v * l for v, l in zip(value_list, (l1, l0))]
    for x, y in zip(*value_list):
        position = dict((k, z) for k, z in zip(posindexes.keys(), [x, y]))
        label = "{:.2f}{} {:.2f}{}".format(lat.isel(position).data, unitlat,
                                   lon.isel(position).data, unitlon)  #fixme: depth if different
        ax.plot(var.isel(position).data, label=label, **kwargs)

    plt.legend()
    plt.tight_layout()

    return ax


if __name__ == '__main__':
    """Tests
    """
    import xarray as xr

    ds = xr.open_dataset('http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_v2/gbr4_simple_2019-01.nc')

    proj = ccrs.PlateCarree()
    fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw=dict(projection=proj))

    plot_map(ds, 'temp', {'time': 0, 'zc': -0.5}, ax=ax1, cmap='inferno')
    ax2 = plot_map(ds, 'salt', {'time': 0, 'k': -3,
                                'j': slice(10, 100), 'i': slice(100)},
                   ax=ax2, vmin=32, add_coastline=False)
    ax2.coastlines(resolution='10m')

    fig.suptitle('Plotting on different scales')

    plt.show()

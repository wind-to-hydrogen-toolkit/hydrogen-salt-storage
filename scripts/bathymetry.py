#!/usr/bin/env python
# coding: utf-8

# # Bathymetry
#
# Elevation relative to sea level /
# Sea-floor height (above Lowest Astronomical Tide datum)
#
# - <https://emodnet.ec.europa.eu/en/bathymetry>
# - <https://doi.org/10.12770/ff3aff8a-cff1-44a3-a2c8-1910bf109f85>
# - <https://emodnet.ec.europa.eu/geonetwork/emodnet/eng/catalog.search#/metadata/53e69177-16cc-4b7a-a6e1-2a4f245e4dbd>

# In[ ]:


import os
from zipfile import BadZipFile, ZipFile

import cartopy.crs as ccrs
import contextily as cx
import matplotlib.pyplot as plt
import numpy as np
import rasterio as rio
import seaborn as sns
import xarray as xr
from matplotlib_scalebar.scalebar import ScaleBar

from h2ss import data as rd

# In[ ]:


plt.rcParams["xtick.major.size"] = 0
plt.rcParams["ytick.major.size"] = 0


# In[ ]:


# base data download directory
DATA_DIR = os.path.join("data", "bathymetry")

# DTM tile D4
FILE_NAME = "D4_2022.nc.zip"

URL = f"https://downloads.emodnet-bathymetry.eu/v11/{FILE_NAME}"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))


# In[ ]:


rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)


# In[ ]:


ZipFile(DATA_FILE).namelist()


# In[ ]:


# extract the archive
try:
    z = ZipFile(DATA_FILE)
    z.extractall(DATA_DIR)
except BadZipFile:
    print("There were issues with the file", DATA_FILE)


# In[ ]:


data = xr.open_dataset(
    os.path.join(DATA_DIR, "D4_2022.nc"), decode_coords="all"
)


# In[ ]:


data = data.chunk({"lat": 1000, "lon": 1000, "cdi_index_count": 1000})


# In[ ]:


data


# In[ ]:


data.rio.crs


# In[ ]:


data.rio.bounds()


# In[ ]:


data.rio.resolution()


# In[ ]:


# read Kish Basin data and extent
ds, extent = rd.read_dat_file(dat_path=os.path.join("data", "kish-basin"))


# In[ ]:


xmin, ymin, xmax, ymax = extent.total_bounds


# In[ ]:


shape = rd.halite_shape(dat_xr=ds).buffer(1000).buffer(-1000)


# In[ ]:


# reproject bathymetry data to the Kish Basin data's CRS
data_ = data.rio.reproject(rd.CRS).rio.clip(extent)


# In[ ]:


data_.rio.resolution()


# In[ ]:


def plot_bath_map(xds, levels=15, cmap="mako", vmax=None, vmin=None):
    """Plotting helper function"""
    plt.figure(figsize=(10, 7))
    ax = plt.axes(projection=ccrs.epsg(rd.CRS))
    xds.plot.contourf(
        ax=ax,
        robust=True,
        cmap=cmap,
        levels=levels,
        vmax=vmax,
        vmin=vmin,
        extend="both",
    )
    shape.boundary.plot(ax=ax, color="white")

    plt.ylim(ymin - 3000, ymax + 3000)
    plt.xlim(xmin - 3000, xmax + 3000)

    cx.add_basemap(
        ax, source=cx.providers.CartoDB.Positron, crs=rd.CRS, zoom=10
    )
    ax.gridlines(
        draw_labels={"bottom": "x", "left": "y"},
        alpha=0.25,
        color="darkslategrey",
    )
    ax.add_artist(
        ScaleBar(1, box_alpha=0, location="lower right", color="darkslategrey")
    )

    plt.title(None)
    plt.tight_layout()
    plt.show()


# In[ ]:


plot_bath_map(data_["elevation"])


# ## Reproject bathymetry to match the resolution of the Kish Basin data

# In[ ]:


data_ = data.rename({"lon": "x", "lat": "y"}).rio.reproject_match(
    ds, resampling=rio.enums.Resampling.bilinear
)


# In[ ]:


data_


# In[ ]:


data_.rio.crs


# In[ ]:


data_.rio.bounds()


# In[ ]:


data_.rio.resolution()


# In[ ]:


val = data_["elevation"].values.flatten()
min(val[~np.isnan(val)])


# In[ ]:


max(val[~np.isnan(val)])


# In[ ]:


plot_bath_map(data_["elevation"], levels=14, vmax=0, vmin=-130)


# ## Adjust Kish Basin depth from sea level to seabed

# In[ ]:


ds = ds.assign(TopDepthSeabed=ds["TopDepth"] + data_["elevation"])
ds = ds.assign(BaseDepthSeabed=ds["BaseDepth"] + data_["elevation"])


# In[ ]:


# ds["TopDepthSeabed"].attrs = ds["TopDepth"].attrs
# ds["BaseDepthSeabed"].attrs = ds["BaseDepth"].attrs


# In[ ]:


ds


# In[ ]:


val = ds["TopDepthSeabed"].values.flatten()
min(val[~np.isnan(val)])


# In[ ]:


max(val[~np.isnan(val)])


# In[ ]:


val = ds["BaseDepthSeabed"].values.flatten()
min(val[~np.isnan(val)])


# In[ ]:


max(val[~np.isnan(val)])


# In[ ]:


plot_bath_map(
    ds["TopDepthSeabed"].sel(halite="Rossall"), cmap="jet", vmin=150, vmax=2250
)


# In[ ]:


plot_bath_map(ds["BaseDepthSeabed"].sel(halite="Rossall"), cmap="jet")


# In[ ]:


bath = rd.bathymetry_layer(
    dat_extent=extent.buffer(3000),
    bathymetry_path=os.path.join("data", "bathymetry"),
)


# In[ ]:


plt.figure(figsize=(10, 7))
ax2 = plt.axes(projection=ccrs.epsg(rd.CRS))
bath["elevation"].plot.contourf(cmap="mako", levels=10, robust=True)
CS = bath["elevation"].plot.contour(
    # colors="darkslategrey",
    colors="white",
    levels=10,
    linewidths=0.5,
    linestyles="solid",
    # alpha=0.5,
    robust=True,
)
plt.clabel(CS, inline=True)
cx.add_basemap(ax2, source=cx.providers.CartoDB.VoyagerNoLabels, crs=rd.CRS)
cx.add_basemap(ax2, source=cx.providers.CartoDB.VoyagerOnlyLabels, crs=rd.CRS)
ax2.gridlines(
    draw_labels={"bottom": "x", "left": "y"},
    alpha=0.25,
    color="darkslategrey",
)
plt.title(None)
plt.tight_layout()
plt.show()

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

# In[1]:


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

# In[2]:


plt.rcParams["xtick.major.size"] = 0
plt.rcParams["ytick.major.size"] = 0


# In[3]:


# base data download directory
DATA_DIR = os.path.join("data", "bathymetry")

# DTM tile D4
FILE_NAME = "D4_2022.nc.zip"

URL = f"https://downloads.emodnet-bathymetry.eu/v11/{FILE_NAME}"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))


# In[4]:


rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)


# In[5]:


ZipFile(DATA_FILE).namelist()


# In[6]:


# extract the archive
try:
    z = ZipFile(DATA_FILE)
    z.extractall(DATA_DIR)
except BadZipFile:
    print("There were issues with the file", DATA_FILE)


# In[7]:


data = xr.open_dataset(
    os.path.join(DATA_DIR, "D4_2022.nc"), decode_coords="all"
)


# In[8]:


data = data.chunk({"lat": 1000, "lon": 1000, "cdi_index_count": 1000})


# In[9]:


data


# In[10]:


data.rio.crs


# In[11]:


data.rio.bounds()


# In[12]:


data.rio.resolution()


# In[4]:


# read Kish Basin data and extent
ds, extent = rd.read_dat_file(dat_path=os.path.join("data", "kish-basin"))


# In[14]:


xmin, ymin, xmax, ymax = extent.total_bounds


# In[15]:


shape = rd.halite_shape(dat_xr=ds).buffer(1000).buffer(-1000)


# In[16]:


# reproject bathymetry data to the Kish Basin data's CRS
data_ = data.rio.reproject(rd.CRS).rio.clip(extent)


# In[17]:


data_.rio.resolution()


# In[18]:


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


# In[19]:


plot_bath_map(data_["elevation"])


# ## Reproject bathymetry to match the resolution of the Kish Basin data

# In[20]:


data_ = data.rename({"lon": "x", "lat": "y"}).rio.reproject_match(
    ds, resampling=rio.enums.Resampling.bilinear
)


# In[21]:


data_


# In[22]:


data_.rio.crs


# In[23]:


data_.rio.bounds()


# In[24]:


data_.rio.resolution()


# In[38]:


val = data_["elevation"].values.flatten()
min(val[~np.isnan(val)])


# In[39]:


max(val[~np.isnan(val)])


# In[27]:


plot_bath_map(data_["elevation"], levels=14, vmax=0, vmin=-130)


# ## Adjust Kish Basin depth from sea level to seabed

# In[28]:


ds = ds.assign(TopDepthSeabed=ds["TopDepth"] + data_["elevation"])
ds = ds.assign(BaseDepthSeabed=ds["BaseDepth"] + data_["elevation"])


# In[29]:


# ds["TopDepthSeabed"].attrs = ds["TopDepth"].attrs
# ds["BaseDepthSeabed"].attrs = ds["BaseDepth"].attrs


# In[30]:


ds


# In[40]:


val = ds["TopDepthSeabed"].values.flatten()
min(val[~np.isnan(val)])


# In[41]:


max(val[~np.isnan(val)])


# In[42]:


val = ds["BaseDepthSeabed"].values.flatten()
min(val[~np.isnan(val)])


# In[43]:


max(val[~np.isnan(val)])


# In[35]:


plot_bath_map(
    ds["TopDepthSeabed"].sel(halite="Rossall"), cmap="jet", vmin=150, vmax=2250
)


# In[36]:


plot_bath_map(ds["BaseDepthSeabed"].sel(halite="Rossall"), cmap="jet")


# In[5]:


bath = rd.bathymetry_layer(
    dat_extent=extent.buffer(3000),
    bathymetry_path=os.path.join("data", "bathymetry"),
)


# In[11]:


plt.figure(figsize=(10, 7))
ax = plt.axes(projection=ccrs.epsg(rd.CRS))
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
cx.add_basemap(ax, source=cx.providers.CartoDB.VoyagerNoLabels, crs=rd.CRS)
cx.add_basemap(ax, source=cx.providers.CartoDB.VoyagerOnlyLabels, crs=rd.CRS)
ax.gridlines(
    draw_labels={"bottom": "x", "left": "y"},
    alpha=0.25,
    color="darkslategrey",
)
plt.title(None)
plt.tight_layout()
plt.show()

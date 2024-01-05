#!/usr/bin/env python
# coding: utf-8

# # Bathymetry
#
# <https://emodnet.ec.europa.eu/en/bathymetry>

import os
from zipfile import BadZipFile, ZipFile

import cartopy.crs as ccrs
import contextily as cx
import matplotlib.pyplot as plt
import rasterio as rio
import rioxarray as rxr
import seaborn as sns
from matplotlib_scalebar.scalebar import ScaleBar

from src import data as rd

# base data download directory
DATA_DIR = os.path.join("data", "bathymetry")

# DTM tile D4
FILE_NAME = "D4_2022.tif.zip"

URL = "https://downloads.emodnet-bathymetry.eu/v11/" + FILE_NAME

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

# extract the archive
try:
    z = ZipFile(DATA_FILE)
    z.extractall(DATA_DIR)
except BadZipFile:
    print("There were issues with the file", DATA_FILE)

data = rxr.open_rasterio(
    os.path.join(DATA_DIR, "D4_2022_mean.tif"), chunks="auto"
)

data

data.rio.crs

data.rio.bounds()

data.rio.resolution()

# read Kish Basin data and extent
ds, extent = rd.read_dat_file(dat_path=os.path.join("data", "kish-basin"))

xmin, ymin, xmax, ymax = extent.total_bounds

shape = rd.halite_shape(dat_xr=ds).buffer(1000).buffer(-1000)

data.rio.reproject(rd.CRS).rio.resolution()

plt.figure(figsize=(10, 7))
ax = plt.axes(projection=ccrs.epsg(rd.CRS))
data.rio.reproject(rd.CRS).rio.clip(extent).isel(band=0).plot.contourf(
    ax=ax, robust=True, cmap="mako", levels=15
)
shape.boundary.plot(ax=ax, color="white")

plt.ylim(ymin - 3000, ymax + 3000)
plt.xlim(xmin - 3000, xmax + 3000)

cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, crs=rd.CRS, zoom=10)
ax.gridlines(
    draw_labels={"bottom": "x", "left": "y"}, alpha=0.25, color="darkslategrey"
)
ax.add_artist(
    ScaleBar(1, box_alpha=0, location="lower right", color="darkslategrey")
)

plt.title(None)
plt.tight_layout()
plt.show()

# ## Reproject bathymetry to match the Kish Basin data

data_ = data.rio.reproject_match(ds, resampling=rio.enums.Resampling.bilinear)

data_

data_.rio.crs

data_.rio.bounds()

data_.rio.resolution()

plt.figure(figsize=(10, 7))
ax = plt.axes(projection=ccrs.epsg(rd.CRS))
data_.isel(band=0).plot.contourf(
    ax=ax, extend="both", cmap="mako", levels=14, vmax=0, vmin=-130
)
shape.boundary.plot(ax=ax, color="white")

plt.ylim(ymin - 3000, ymax + 3000)
plt.xlim(xmin - 3000, xmax + 3000)

cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, crs=rd.CRS, zoom=10)
ax.gridlines(
    draw_labels={"bottom": "x", "left": "y"}, alpha=0.25, color="darkslategrey"
)
ax.add_artist(
    ScaleBar(1, box_alpha=0, location="lower right", color="darkslategrey")
)

plt.title(None)
plt.tight_layout()
plt.show()

# ## Adjust Kish Basin depth from sea level to seabed

ds = ds.assign(TopDepthSeabed=ds["TopDepth"] + data_.isel(band=0))
ds = ds.assign(BaseDepthSeabed=ds["BaseDepth"] + data_.isel(band=0))

ds["TopDepthSeabed"].attrs = ds["TopDepth"].attrs
ds["BaseDepthSeabed"].attrs = ds["BaseDepth"].attrs

ds

plt.figure(figsize=(10, 9))
ax = plt.axes(projection=ccrs.epsg(rd.CRS))
data_.isel(band=0).plot.contourf(
    ax=ax,
    extend="both",
    cmap="mako",
    levels=14,
    vmax=0,
    vmin=-130,
    cbar_kwargs={"pad": 0.045, "shrink": 0.8, "aspect": 25},
)
shape.boundary.plot(ax=ax, color="white")

ds["TopDepthSeabed"].sel(halite="Rossall").plot.contourf(
    cmap="jet",
    levels=15,
    robust=True,
    ax=ax,
    cbar_kwargs={"location": "bottom", "aspect": 30, "pad": 0.045},
)

plt.ylim(ymin - 3000, ymax + 3000)
plt.xlim(xmin - 3000, xmax + 3000)

cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, crs=rd.CRS, zoom=10)
ax.gridlines(
    draw_labels={"bottom": "x", "left": "y"}, alpha=0.25, color="darkslategrey"
)
ax.add_artist(
    ScaleBar(1, box_alpha=0, location="lower right", color="darkslategrey")
)

plt.title(None)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 9))
ax = plt.axes(projection=ccrs.epsg(rd.CRS))
data_.isel(band=0).plot.contourf(
    ax=ax,
    extend="both",
    cmap="mako",
    levels=14,
    vmax=0,
    vmin=-130,
    cbar_kwargs={"pad": 0.045, "shrink": 0.8, "aspect": 25},
)
shape.boundary.plot(ax=ax, color="white")

ds["BaseDepthSeabed"].sel(halite="Rossall").plot.contourf(
    cmap="jet",
    levels=15,
    robust=True,
    ax=ax,
    cbar_kwargs={"location": "bottom", "aspect": 30, "pad": 0.045},
)

plt.ylim(ymin - 3000, ymax + 3000)
plt.xlim(xmin - 3000, xmax + 3000)

cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, crs=rd.CRS, zoom=10)
ax.gridlines(
    draw_labels={"bottom": "x", "left": "y"}, alpha=0.25, color="darkslategrey"
)
ax.add_artist(
    ScaleBar(1, box_alpha=0, location="lower right", color="darkslategrey")
)

plt.title(None)
plt.tight_layout()
plt.show()

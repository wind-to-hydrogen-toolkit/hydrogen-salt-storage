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

import os
from zipfile import BadZipFile, ZipFile

import cartopy.crs as ccrs
import contextily as cx
import matplotlib.pyplot as plt
import rasterio as rio
import seaborn as sns
import xarray as xr
from matplotlib_scalebar.scalebar import ScaleBar

from h2ss import data as rd

# base data download directory
DATA_DIR = os.path.join("data", "bathymetry")

# DTM tile D4
FILE_NAME = "D4_2022.nc.zip"

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

data = xr.open_dataset(
    os.path.join(DATA_DIR, "D4_2022.nc"), decode_coords="all"
)

data = data.chunk({"lat": 1000, "lon": 1000, "cdi_index_count": 1000})

data

data.rio.crs

data.rio.bounds()

data.rio.resolution()

# read Kish Basin data and extent
ds, extent = rd.read_dat_file(dat_path=os.path.join("data", "kish-basin"))

xmin, ymin, xmax, ymax = extent.total_bounds

shape = rd.halite_shape(dat_xr=ds).buffer(1000).buffer(-1000)

# reproject bathymetry data to the Kish Basin data's CRS
data_ = data.rio.reproject(rd.CRS).rio.clip(extent)

data_.rio.resolution()


def plot_bath_map(xds, levels=15, cmap="mako", vmax=None, vmin=None):
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


plot_bath_map(data_["elevation"])

# ## Reproject bathymetry to match the resolution of the Kish Basin data

data_ = data.rename({"lon": "x", "lat": "y"}).rio.reproject_match(
    ds, resampling=rio.enums.Resampling.bilinear
)

data_

data_.rio.crs

data_.rio.bounds()

data_.rio.resolution()

min(set(data_["elevation"].values.flatten()))

max(set(data_["elevation"].values.flatten()))

plot_bath_map(data_["elevation"], levels=14, vmax=0, vmin=-130)

# ## Adjust Kish Basin depth from sea level to seabed

ds = ds.assign(TopDepthSeabed=ds["TopDepth"] + data_["elevation"])
ds = ds.assign(BaseDepthSeabed=ds["BaseDepth"] + data_["elevation"])

ds["TopDepthSeabed"].attrs = ds["TopDepth"].attrs
ds["BaseDepthSeabed"].attrs = ds["BaseDepth"].attrs

ds

plot_bath_map(ds["TopDepthSeabed"].sel(halite="Rossall"), cmap="jet")

plot_bath_map(ds["BaseDepthSeabed"].sel(halite="Rossall"), cmap="jet")

min(set(ds["TopDepthSeabed"].values.flatten()))

max(set(ds["TopDepthSeabed"].values.flatten()))

min(set(ds["BaseDepthSeabed"].values.flatten()))

max(set(ds["BaseDepthSeabed"].values.flatten()))

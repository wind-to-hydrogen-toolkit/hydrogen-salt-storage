#!/usr/bin/env python
# coding: utf-8

# # Kish Basin Salt Caverns
#
# <https://hyss.ie/>

import os
from zipfile import BadZipFile, ZipFile
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import pooch
import glob
from datetime import datetime, timezone
from shapely.geometry import Polygon
import rioxarray as rxr
from geocube.api.core import make_geocube
from textwrap import wrap
import numpy as np
import shapely
import itertools

# base data download directory
DATA_DIR = os.path.join("data", "kish-basin")
FILE_NAME = "Kish-Basin-dat-files.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# boundary data
ie = gpd.read_file(
    os.path.join("data", "boundaries.gpkg"), layer="NUTS_RG_01M_2021_4326_IE"
)

crs = 23029

# ## Read data layers


def read_dat_file(dat_path: str, dat_crs):
    """
    Read XYZ data layers into an Xarray dataset
    """

    gdf = {}
    for dat_file in glob.glob(dat_path):
        # read each layer as individual dataframes into a dictionary
        gdf[os.path.split(dat_file)[1][:-4]] = pd.read_fwf(
            dat_file, header=None, names=["X", "Y", "Z"]
        )

        # assign layer name to a column
        gdf[os.path.split(dat_file)[1][:-4]]["data"] = os.path.split(dat_file)[
            1
        ][:-4]

    # combine dataframes
    gdf = pd.concat(gdf.values())

    # convert dataframe to geodataframe
    gdf["wkt"] = (
        "POINT (" + gdf["X"].astype(str) + " " + gdf["Y"].astype(str) + ")"
    )
    gdf = gpd.GeoDataFrame(
        gdf, geometry=gpd.GeoSeries.from_wkt(gdf["wkt"]), crs=dat_crs
    )
    gdf.drop(columns=["wkt", "X", "Y"], inplace=True)

    # convert to Xarray dataset
    ds = make_geocube(vector_data=gdf, resolution=(200, -200), group_by="data")

    # keep only zones of interest in the dataframe
    gdf = gdf.loc[gdf["data"].str.contains("Zone")]

    return gdf, ds


gdf, ds = read_dat_file(os.path.join(DATA_DIR, "*.dat"), dat_crs=crs)

# ### Map extent

# create extent polygon
extent = pd.read_csv(
    os.path.join(DATA_DIR, "Kish GIS Map Extent - Square.csv"), skiprows=2
)
extent = gpd.GeoSeries(
    Polygon(
        [
            (extent[" X"][0], extent[" Y"][0]),
            (extent[" X"][1], extent[" Y"][1]),
            (extent[" X"][2], extent[" Y"][2]),
            (extent[" X"][3], extent[" Y"][3]),
        ]
    ),
    crs=crs,
)

extent.bounds

# ### Halite layers

ds

ds.rio.crs

ds.rio.resolution()

ds.rio.bounds()

# ### Zone of interest boundaries

zones = ds.sel(data=[x for x in ds["data"].values if "Zone" in x])

zones

zones.rio.bounds()

xmin, ymin, xmax, ymax = zones.rio.bounds()

cell_size = 200
# create the cells in a loop
grid_cells = []
for x0 in np.arange(xmin, xmax + cell_size, cell_size):
    for y0 in np.arange(ymin, ymax + cell_size, cell_size):
        # bounds
        x1 = x0 - cell_size
        y1 = y0 + cell_size
        grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))
grid_cells = gpd.GeoDataFrame(grid_cells, columns=["geometry"], crs=crs)

# drop cells without data
grid_centroids = {"wkt": [], "x": [], "y": []}

for x, y in itertools.product(
    range(len(zones.coords["x"])), range(len(zones.coords["y"]))
):
    data__ = zones.isel(x=x, y=y)

    # ignore null cells
    if not data__["Z"].isnull().all():
        grid_centroids["wkt"].append(
            f"POINT ({float(data__['x'].values)} "
            f"{float(data__['y'].values)})"
        )
        grid_centroids["x"].append(float(data__["x"].values))
        grid_centroids["y"].append(float(data__["y"].values))

grid_centroids = gpd.GeoDataFrame(
    grid_centroids,
    geometry=gpd.GeoSeries.from_wkt(grid_centroids["wkt"], crs=crs),
)

zones = gpd.sjoin(grid_cells, grid_centroids)

zones.drop(columns=["wkt", "index_right", "x", "y"], inplace=True)

zones["index_"] = "0"

zones = zones.dissolve(by="index_").reset_index().drop(columns=["index_"])

zones

ax = zones.boundary.plot(linewidth=0.5)
extent.boundary.plot(ax=ax, color="darkslategrey")
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()

# ### Zone of interest points

gdf.shape

gdf.head()

ax = gdf.plot(markersize=0.25)
extent.boundary.plot(ax=ax, color="darkslategrey")
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# ### Zone of interest boundaries

zones = gpd.GeoSeries(gdf["geometry"].buffer(200).unary_union, crs=crs)

zones

zones.bounds

ax = zones.boundary.plot(linewidth=1)
extent.boundary.plot(ax=ax, color="darkslategrey")
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# ### Cavern specifications

# cavern specifications
diameter = 84
separation = diameter * 4

separation

# ### Generate salt cavern grid using extent

# xmin, ymin, xmax, ymax = (
#     extent.bounds["minx"][0], extent.bounds["miny"][0],
#     extent.bounds["maxx"][0], extent.bounds["maxy"][0]
# )

# create the cells in a loop
grid_cells = []
for x0 in np.arange(xmin, xmax + separation, separation):
    for y0 in np.arange(ymin, ymax + separation, separation):
        # bounds
        x1 = x0 - separation
        y1 = y0 + separation
        grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))
grid_cells = gpd.GeoDataFrame(grid_cells, columns=["geometry"], crs=crs)

# verify separation distance
x0 - x1, y0 - y1

grid_cells.head()

grid_cells.shape

ax = grid_cells.boundary.plot(linewidth=0.5)
zones.boundary.plot(ax=ax, linewidth=1, color="darkslategrey")
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# ### Caverns in zone of interest

caverns = gpd.GeoDataFrame(
    geometry=grid_cells.centroid.buffer(diameter / 2)
).overlay(gpd.GeoDataFrame(geometry=zones["geometry"]), how="intersection")

caverns.head()

caverns.shape

ax = caverns.centroid.plot(markersize=1)
zones.boundary.plot(ax=ax, linewidth=1, color="darkslategrey")
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

plt.figure(figsize=(9, 7))
ax = plt.axes(projection=ccrs.epsg(crs))
ds.sel(data="Rossall Halite Thickness - Zone Of Interest - XYZ Meters")[
    "Z"
].plot.contourf(
    ax=ax,
    cmap="jet",
    alpha=0.25,
    robust=True,
    levels=[300 + 30 * n for n in range(15)],
    transform=ccrs.epsg(crs),
    xlim=(687000, 742000),
    ylim=(5888000, 5937000),
    extend="max",
)
ds.sel(
    data="Presall Halite Thickness - Zone Of Interest - XYZ Meters-corrected"
)["Z"].plot.contourf(
    ax=ax,
    cmap="jet",
    alpha=0.25,
    robust=True,
    levels=[300 + 30 * n for n in range(15)],
    transform=ccrs.epsg(crs),
    xlim=(687000, 742000),
    ylim=(5888000, 5937000),
    extend="max",
    add_colorbar=False,
)
ds.sel(data="Flyde Halite Thickness - Zone Of Interest - XYZ Meters")[
    "Z"
].plot.contourf(
    ax=ax,
    cmap="jet",
    alpha=0.25,
    robust=True,
    levels=[300 + 30 * n for n in range(15)],
    transform=ccrs.epsg(crs),
    xlim=(687000, 742000),
    ylim=(5888000, 5937000),
    extend="max",
    add_colorbar=False,
)
caverns.centroid.plot(ax=ax, markersize=2, marker=".", color="black")
ie.to_crs(crs).boundary.plot(ax=ax, edgecolor="darkslategrey", linewidth=0.5)
plt.title("")
plt.tight_layout()
plt.show()

plt.figure(figsize=(9, 7))
ax = plt.axes(projection=ccrs.epsg(crs))
ds.sel(data="Rossall Halite Thickness XYZ Meters")["Z"].plot.contourf(
    ax=ax,
    cmap="jet",
    alpha=0.25,
    robust=True,
    levels=15,
    transform=ccrs.epsg(crs),
    xlim=(687000, 742000),
    ylim=(5888000, 5937000),
)
caverns.centroid.plot(ax=ax, markersize=2, marker=".", color="black")
ie.to_crs(crs).boundary.plot(ax=ax, edgecolor="darkslategrey", linewidth=0.5)
plt.title("")
plt.tight_layout()
plt.show()

#!/usr/bin/env python
# coding: utf-8

# # Kish Basin Salt Caverns
#
# <https://hyss.ie/>

import glob
import itertools
import os

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray as rxr
import shapely
from geocube.api.core import make_geocube

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
    gdf = make_geocube(
        vector_data=gdf, resolution=(200, -200), group_by="data"
    )

    return gdf


ds = read_dat_file(os.path.join(DATA_DIR, "*.dat"), dat_crs=crs)

# ### Map extent

# create extent polygon
extent = pd.read_csv(
    os.path.join(DATA_DIR, "Kish GIS Map Extent - Square.csv"), skiprows=2
)
extent = gpd.GeoSeries(
    shapely.geometry.Polygon(
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

# ### Zones of interest boundaries

zones = ds.sel(data=[x for x in ds["data"].values if "Zone" in x])

zones

zones.rio.bounds()

# create gridded layer
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
    data = zones.isel(x=x, y=y)

    # ignore null cells
    if not data["Z"].isnull().all():
        grid_centroids["wkt"].append(
            f"POINT ({float(data['x'].values)} " f"{float(data['y'].values)})"
        )
        grid_centroids["x"].append(float(data["x"].values))
        grid_centroids["y"].append(float(data["y"].values))

grid_centroids = gpd.GeoDataFrame(
    grid_centroids,
    geometry=gpd.GeoSeries.from_wkt(grid_centroids["wkt"], crs=crs),
)

zones = gpd.sjoin(grid_cells, grid_centroids)
zones.drop(columns=["wkt", "index_right", "x", "y"], inplace=True)
zones["index_"] = "0"

# convert gridded data to polygon
zones = zones.dissolve(by="index_").reset_index().drop(columns=["index_"])

zones

ax = zones.plot(linewidth=0.5)
extent.boundary.plot(ax=ax, color="darkslategrey")
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()

# ### Cavern specifications

# cavern specifications
diameter = 84
separation = diameter * 4

separation

# ### Generate salt cavern grid using extent

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

ax = grid_cells.boundary.plot(linewidth=0.5, color="darkslategrey")
zones.plot(ax=ax, linewidth=1)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# ### Caverns in zones of interest

caverns = gpd.GeoDataFrame(
    geometry=grid_cells.centroid.buffer(diameter / 2)
).overlay(gpd.GeoDataFrame(geometry=zones["geometry"]), how="intersection")

caverns.head()

caverns.shape[0]

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
    cmap="jet",
    alpha=0.25,
    levels=[300 + 30 * n for n in range(15)],
    xlim=(extent.bounds["minx"][0], extent.bounds["maxx"][0]),
    ylim=(extent.bounds["miny"][0], extent.bounds["maxy"][0]),
    extend="max",
    cbar_kwargs={"label": "Halite Thickness (m)"},
)
ds.sel(
    data="Presall Halite Thickness - Zone Of Interest - XYZ Meters-corrected"
)["Z"].plot.contourf(
    cmap="jet",
    alpha=0.25,
    levels=[300 + 30 * n for n in range(15)],
    extend="max",
    add_colorbar=False,
)
ds.sel(data="Flyde Halite Thickness - Zone Of Interest - XYZ Meters")[
    "Z"
].plot.contourf(
    cmap="jet",
    alpha=0.25,
    levels=[300 + 30 * n for n in range(15)],
    extend="max",
    add_colorbar=False,
)
caverns.centroid.plot(ax=ax, markersize=2, marker=".", color="black")
ie.to_crs(crs).boundary.plot(ax=ax, edgecolor="darkslategrey", linewidth=0.5)
plt.title("Kish Basin Caverns within Halite Zones of Interest")
plt.tight_layout()
plt.show()

plt.figure(figsize=(9, 7))
ax = plt.axes(projection=ccrs.epsg(crs))
ds.sel(data="Rossall Halite Thickness XYZ Meters")["Z"].plot.contourf(
    cmap="jet",
    alpha=0.25,
    robust=True,
    levels=15,
    xlim=(extent.bounds["minx"][0], extent.bounds["maxx"][0]),
    ylim=(extent.bounds["miny"][0], extent.bounds["maxy"][0]),
    cbar_kwargs={"label": "Halite Thickness (m)"},
)
caverns.centroid.plot(ax=ax, markersize=2, marker=".", color="black")
ie.to_crs(crs).boundary.plot(ax=ax, edgecolor="darkslategrey", linewidth=0.5)
plt.title("Kish Basin Caverns within Rossall Halite")
plt.tight_layout()
plt.show()

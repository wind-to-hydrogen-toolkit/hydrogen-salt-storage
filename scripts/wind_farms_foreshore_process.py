#!/usr/bin/env python
# coding: utf-8

# # Wind Farms (Foreshore Process)
#
# <https://data.gov.ie/dataset/wind-farms-foreshore-process>

import os
from datetime import datetime, timezone
from zipfile import ZipFile

import cartopy.crs as ccrs
import contextily as cx
import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pooch
from matplotlib_scalebar.scalebar import ScaleBar

from src import functions as fns

# base data download directory
DATA_DIR = os.path.join("data", "wind-farms-foreshore-process")
os.makedirs(DATA_DIR, exist_ok=True)

URL = (
    "https://opendata.arcgis.com/api/v3/datasets/"
    "803a4ecc22aa4cc09111072a0bbc4fac_2/downloads/"
    "data?format=shp&spatialRefId=4326&where=1%3D1"
)
KNOWN_HASH = None
FILE_NAME = "wind-farms-foreshore-process.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

# download data if necessary
if not os.path.isfile(DATA_FILE):
    pooch.retrieve(
        url=URL, known_hash=KNOWN_HASH, fname=FILE_NAME, path=DATA_DIR
    )

    with open(f"{DATA_FILE[:-4]}.txt", "w", encoding="utf-8") as outfile:
        outfile.write(
            f"Data downloaded on: {datetime.now(tz=timezone.utc)}\n"
            f"Download URL: {URL}"
        )

with open(f"{DATA_FILE[:-4]}.txt") as f:
    print(f.read())

ZipFile(DATA_FILE).namelist()

wind_farms = gpd.read_file(
    os.path.join(
        f"zip://{DATA_FILE}!"
        + [x for x in ZipFile(DATA_FILE).namelist() if x.endswith(".shp")][0]
    )
)

wind_farms.crs

wind_farms.shape

wind_farms.columns

wind_farms[["Name", "Type", "MDM_Catego"]]

# minor name fix
wind_farms.at[1, "Name"] = "Kilmichael Point"

ax = wind_farms.to_crs(3857).plot(
    column="Name",
    cmap="tab20",
    alpha=0.5,
    figsize=(10, 10),
    legend=True,
    legend_kwds={"loc": "upper right"},
    linewidth=0.5,
    edgecolor="darkslategrey",
)
plt.xlim(-1.2e6, -0.3e6)
plt.ylim(6.65e6, 7.475e6)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=7)

plt.title("Wind Farms (Foreshore Process)")

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# read Kish Basin data
DATA_DIR = os.path.join("data", "kish-basin")

CRS = 23029

ds, extent = fns.read_dat_file(DATA_DIR, CRS)

# use extent bounds
xmin, ymin, xmax, ymax = extent.total_bounds

# wind farms in the Irish Sea
wind_farms_ = wind_farms.sjoin(
    gpd.GeoDataFrame(geometry=extent.buffer(50000)).to_crs(wind_farms.crs)
)

ax = wind_farms_.to_crs(3857).plot(
    column="Name",
    alpha=0.5,
    figsize=(9, 9),
    cmap="tab20",
    linewidth=0.5,
    edgecolor="darkslategrey",
    legend=True,
)
plt.xlim(-7.5e5, -5.25e5)
extent.to_crs(3857).boundary.plot(ax=ax, alpha=0)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)

plt.title("Wind Farms in the Irish Sea")
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# wind farms near Kish Basin
wind_farms_ = (
    wind_farms.sjoin(
        gpd.GeoDataFrame(geometry=extent.buffer(3000)).to_crs(wind_farms.crs)
    )
    .reset_index()
    .sort_values("Name")
)

plt.figure(figsize=(12, 9))
ax = plt.axes(projection=ccrs.epsg(CRS))

ds.max(dim="halite")["Thickness"].plot.contourf(
    cmap="jet",
    alpha=0.65,
    robust=True,
    levels=15,
    cbar_kwargs={"label": "Maximum Halite Thickness [m]"},
)

plt.xlim(xmin, xmax)
plt.ylim(ymin - 4000, ymax + 4000)

# wind farms
colours = ["firebrick", "forestgreen", "black", "royalblue"]
for index, colour in zip(range(len(wind_farms_)), colours):
    wind_farms_.iloc[[index]].to_crs(CRS).plot(
        ax=ax, hatch="///", facecolor="none", edgecolor=colour, linewidth=2
    )
legend_handles = [
    mpatches.Patch(
        facecolor="none",
        hatch="////",
        edgecolor=colours[x],
        label=list(wind_farms_["Name"])[x],
    )
    for x in range(len(wind_farms_))
]

cx.add_basemap(ax, crs=CRS, source=cx.providers.CartoDB.Voyager)
ax.gridlines(
    draw_labels={"bottom": "x", "left": "y"}, alpha=0.25, color="darkslategrey"
)
ax.add_artist(
    ScaleBar(1, box_alpha=0, location="lower right", color="darkslategrey")
)
ax.legend(handles=legend_handles, loc="lower right", bbox_to_anchor=(1, 0.05))

# plt.title("Wind Farms near Kish Basin")
plt.title(None)
plt.tight_layout()
plt.show()

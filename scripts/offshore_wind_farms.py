#!/usr/bin/env python
# coding: utf-8

# # Offshore wind farms and pipelines
#
# - <https://data.gov.ie/dataset/marine-area-consent-wind>
# - <https://data-housinggovie.opendata.arcgis.com/datasets/housinggovie::marine-area-consent-wind>
# - <https://data.gov.ie/dataset/offshore-gas-pipelines>
# - <https://data-housinggovie.opendata.arcgis.com/datasets/housinggovie::offshore-gas-pipelines>

import os
from zipfile import ZipFile

import cartopy.crs as ccrs
import contextily as cx
import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from matplotlib_scalebar.scalebar import ScaleBar

from h2ss import data as rd

# base data download directory
DATA_DIR = os.path.join("data", "wind-farms")

URL = (
    "https://opendata.arcgis.com/api/v3/datasets/"
    "803a4ecc22aa4cc09111072a0bbc4fac_3/downloads/"
    "data?format=shp&spatialRefId=4326&where=1%3D1"
)

FILE_NAME = "marine-area-consent-wind.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

plt.rcParams["xtick.major.size"] = 0
plt.rcParams["ytick.major.size"] = 0
plt.rcParams["xtick.minor.size"] = 0
plt.rcParams["ytick.minor.size"] = 0

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

wind_farms = rd.read_shapefile_from_zip(data_path=DATA_FILE)

wind_farms.crs

wind_farms.shape

wind_farms.columns

wind_farms

ax = wind_farms.to_crs(3857).plot(
    column="name",
    cmap="tab20",
    alpha=0.5,
    figsize=(10, 10),
    legend=True,
    legend_kwds={"loc": "lower right"},
    linewidth=0.5,
    edgecolor="darkslategrey",
)
plt.xlim(-1.2e6, -0.55e6)
plt.ylim(6.65e6, 7.475e6)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=7)

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# keep only wind farms near the Kish Basin
wind_farms.drop(index=[0, 1, 7], inplace=True)

# merge wind farm polygons
wind_farms.at[2, "name"] = "Dublin Array"
wind_farms.at[3, "name"] = "Dublin Array"
wind_farms = wind_farms.dissolve(by="name")
wind_farms.reset_index(inplace=True)

wind_farms

ax = wind_farms.to_crs(3857).plot(
    column="name",
    alpha=0.5,
    figsize=(10, 10),
    legend=True,
    legend_kwds={"loc": "upper right"},
    linewidth=0.5,
    edgecolor="darkslategrey",
)
plt.xlim(-0.72e6, -0.62e6)
plt.ylim(6.975e6, 7.15e6)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# read Kish Basin data
DATA_DIR = os.path.join("data", "kish-basin")

ds, extent = rd.read_dat_file(dat_path=DATA_DIR)

# use extent bounds
xmin, ymin, xmax, ymax = extent.total_bounds

# shape of the halite
shape = rd.halite_shape(dat_xr=ds)

plt.figure(figsize=(11, 11))
ax = plt.axes(projection=ccrs.epsg(rd.CRS))

ds.max(dim="halite")["Thickness"].plot.contourf(
    cmap="jet",
    alpha=0.65,
    robust=True,
    levels=15,
    cbar_kwargs={
        "label": "Maximum halite thickness [m]",
        "aspect": 25,
        "pad": 0.035,
    },
)

plt.xlim(xmin - 8550, xmax + 1000)
# plt.ylim(ymin - 10500, ymax + 10500)

# wind farms
colours = ["firebrick", "black", "royalblue"]
legend_handles = []
for index, colour in zip(range(len(wind_farms)), colours):
    wind_farms.iloc[[index]].to_crs(rd.CRS).to_crs(rd.CRS).plot(
        ax=ax, hatch="///", facecolor="none", linewidth=2, edgecolor=colour
    )
    legend_handles.append(
        mpatches.Patch(
            facecolor="none",
            hatch="///",
            edgecolor=colour,
            label=wind_farms.iloc[[index]]["name"].values[0],
        )
    )

cx.add_basemap(ax, crs=rd.CRS, source=cx.providers.CartoDB.Voyager)
ax.gridlines(
    draw_labels={"bottom": "x", "left": "y"},
    alpha=0.25,
    color="darkslategrey",
    xformatter=LongitudeFormatter(auto_hide=False, dms=True),
    yformatter=LatitudeFormatter(auto_hide=False, dms=True),
    ylabel_style={"rotation": 90},
)
ax.add_artist(
    ScaleBar(
        1,
        box_alpha=0,
        location="lower right",
        color="darkslategrey",
        width_fraction=0.0075,
    )
)
ax.legend(handles=legend_handles, loc="lower right", bbox_to_anchor=(1, 0.05))

plt.title(None)
plt.tight_layout()

# plt.savefig(
#     os.path.join("graphics", "fig_offshore_wind_farms.jpg"),
#     format="jpg",
#     dpi=600,
# )
plt.show()

plt.figure(figsize=(9, 9))
ax = plt.axes(projection=ccrs.epsg(rd.CRS))

# add halite boundary - use buffering to smooth the outline
shape.buffer(1000).buffer(-1000).boundary.plot(
    ax=ax, color="black", linewidth=2
)

plt.xlim(xmin - 7750, xmax + 7750)
plt.ylim(ymin - 10500, ymax + 10500)

# wind farms
colours = ["firebrick", "seagreen", "royalblue"]
legend_handles = []
for index, colour in zip(range(len(wind_farms)), colours):
    wind_farms.iloc[[index]].to_crs(rd.CRS).to_crs(rd.CRS).plot(
        ax=ax, hatch="///", facecolor="none", linewidth=2, edgecolor=colour
    )
    legend_handles.append(
        mpatches.Patch(
            facecolor="none",
            hatch="///",
            edgecolor=colour,
            label=wind_farms.iloc[[index]]["name"].values[0],
        )
    )

legend_handles.append(
    mpatches.Patch(
        facecolor="none",
        edgecolor="black",
        label="Kish Basin halite",
        linewidth=2,
    )
)

cx.add_basemap(ax, crs=rd.CRS, source=cx.providers.CartoDB.Voyager, zoom=10)
ax.gridlines(
    draw_labels={"bottom": "x", "left": "y"},
    alpha=0.25,
    color="darkslategrey",
    xformatter=LongitudeFormatter(auto_hide=False, dms=True),
    yformatter=LatitudeFormatter(auto_hide=False, dms=True),
    ylabel_style={"rotation": 90},
)
ax.add_artist(
    ScaleBar(1, box_alpha=0, location="lower right", color="darkslategrey")
)
ax.legend(handles=legend_handles, loc="lower right", bbox_to_anchor=(1, 0.05))

plt.title(None)
plt.tight_layout()
plt.show()

# distance from Kish Bank
wind_farms_ = wind_farms.to_crs(rd.CRS)

for i in range(len(wind_farms_)):
    pv = (
        wind_farms_.iloc[[i]]
        .distance(shape["geometry"], align=False)
        .values[0]
        / 1000
    )
    print(
        wind_farms_.iloc[[i]]["name"].values[0],
        f"is {pv:,.2f} km away from the Kish Bank",
    )

# ## Distance from pipelines

DATA_DIR = os.path.join("data", "pipelines")

URL = (
    "https://opendata.arcgis.com/api/v3/datasets/"
    "dc6e3849b9fc43bb93078c5d0093bf6a_1/downloads/data?"
    "format=shp&spatialRefId=4326&where=1%3D1"
)

FILE_NAME = "pipelines.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

pipelines = rd.read_shapefile_from_zip(data_path=DATA_FILE)

pipelines.head()

pipelines.crs

ax = (
    pipelines.to_crs(rd.CRS)
    .overlay(gpd.GeoDataFrame(geometry=shape.buffer(50000)))
    .plot(color="crimson")
)
shape.buffer(1000).buffer(-1000).boundary.plot(
    ax=ax, color="black", linewidth=1
)
wind_farms_.plot(ax=ax, alpha=0.5, column="name")
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, crs=rd.CRS)

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

for i in range(len(wind_farms_)):
    pv = (
        wind_farms_.iloc[[i]]
        .distance(
            pipelines.to_crs(rd.CRS)
            .overlay(gpd.GeoDataFrame(geometry=shape.buffer(25000)))
            .dissolve()["geometry"],
            align=False,
        )
        .values[0]
        / 1000
    )
    print(
        wind_farms_.iloc[[i]]["name"].values[0],
        f"is {pv:,.2f} km away from the nearest offshore pipeline",
    )

pv = (
    shape.distance(
        pipelines.to_crs(rd.CRS)
        .overlay(gpd.GeoDataFrame(geometry=shape.buffer(25000)))
        .dissolve()["geometry"],
        align=False,
    ).values[0]
    / 1000
)
print(f"Kish Basin is {pv:,.2f} km away from the nearest offshore pipeline")

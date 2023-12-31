#!/usr/bin/env python
# coding: utf-8

# # Weibull Parameters of Wind Speeds 2001 to 2010 -- 150m above ground level
#
# <https://data.gov.ie/dataset/weibull-parameters-wind-speeds-2001-to-2010-150m-above-ground-level>

import os
from datetime import datetime, timezone
from zipfile import ZipFile

import cartopy.crs as ccrs
import contextily as cx
import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pooch
import seaborn as sns
from matplotlib_scalebar.scalebar import ScaleBar

from src import functions as fns

# base data download directory
DATA_DIR = os.path.join("data", "weibull-parameters-wind-speeds")
os.makedirs(DATA_DIR, exist_ok=True)

URL = (
    "https://seaiopendata.blob.core.windows.net/wind/"
    "Weibull_150m_params_ITM.zip"
)
KNOWN_HASH = None
FILE_NAME = "Weibull_150m_params_ITM.zip"

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

with open(f"{DATA_FILE[:-4]}.txt", encoding="utf-8") as f:
    print(f.read())

ZipFile(DATA_FILE).namelist()

weibull_c = fns.read_shapefile_from_zip(
    data_path=os.path.join(DATA_FILE), endswith="c_ITM.shp"
)

weibull_k = fns.read_shapefile_from_zip(
    data_path=os.path.join(DATA_FILE), endswith="k_ITM.shp"
)

weibull_c.crs

weibull_k.crs

weibull_c.shape

weibull_k.shape

weibull_c.columns

weibull_k.columns

weibull_c.head()

weibull_k.head()

ds, extent = fns.read_dat_file(dat_path=os.path.join("data", "kish-basin"))

xmin, ymin, xmax, ymax = extent.total_bounds

ax = weibull_c.to_crs(3857).plot(
    column="Value",
    cmap="flare",
    figsize=(6, 6),
    legend=True,
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

ax = weibull_k.to_crs(3857).plot(
    column="Value",
    cmap="flare",
    figsize=(6, 6),
    legend=True,
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# wind farms in the area
wind_farms = fns.constraint_wind_farm(
    data_path=os.path.join(
        "data",
        "wind-farms-foreshore-process",
        "wind-farms-foreshore-process.zip",
    ),
    dat_extent=extent,
)

# combine Codling wind farm polygons
wind_farms["Name_"] = wind_farms["Name"].str.split(expand=True)[0]

# shape of the halite
shape = fns.halite_shape(dat_xr=ds)

# land boundary
land = fns.read_shapefile_from_zip(
    data_path=os.path.join(
        "data", "boundaries", "osi-provinces-ungeneralised-2019.zip"
    )
)

land = land.dissolve().to_crs(fns.CRS)

# crop to wind farm and basin extent
extent_wf = gpd.GeoDataFrame(geometry=wind_farms.dissolve().envelope)
extent_wf = gpd.GeoDataFrame(
    geometry=(
        extent_wf.overlay(gpd.GeoDataFrame(geometry=extent), how="union")
        .dissolve()
        .envelope
    )
)
weibull_c = weibull_c.to_crs(fns.CRS).overlay(extent_wf, how="intersection")
weibull_k = weibull_k.to_crs(fns.CRS).overlay(extent_wf, how="intersection")

# crop land boundary from c and k
weibull_c = weibull_c.overlay(land, how="difference")
weibull_k = weibull_k.overlay(land, how="difference")

plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.epsg(fns.CRS))

# add halite boundary - use buffering to smooth the outline
shape.buffer(1000).buffer(-1000).boundary.plot(
    ax=ax, color="black", linewidth=2
)

# wind farms
colours = ["lime", "lime", "black", "deepskyblue"]
for index, colour in zip(range(len(wind_farms)), colours):
    wind_farms.iloc[[index]].to_crs(fns.CRS).plot(
        ax=ax,
        hatch="///",
        facecolor="none",
        edgecolor=colour,
        linewidth=2,
        zorder=2,
    )
legend_handles = [
    mpatches.Patch(
        facecolor="none",
        hatch="////",
        edgecolor=colours[1:][x],
        label=list(wind_farms.dissolve("Name_")["Name"])[x],
    )
    for x in range(len(wind_farms.dissolve("Name_")))
]

legend_handles.append(
    mpatches.Patch(
        facecolor="none",
        edgecolor="black",
        label="Kish Basin halite",
        linewidth=2,
    )
)

weibull_c.to_crs(fns.CRS).plot(
    column="Value", cmap="flare", figsize=(6, 6), legend=True, ax=ax, zorder=1
)

cx.add_basemap(ax, crs=fns.CRS, source=cx.providers.CartoDB.Voyager, zoom=10)
ax.gridlines(
    draw_labels={"bottom": "x", "left": "y"}, alpha=0.25, color="darkslategrey"
)
ax.add_artist(
    ScaleBar(1, box_alpha=0, location="lower right", color="darkslategrey")
)
ax.legend(handles=legend_handles, loc="upper right")

plt.title(None)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.epsg(fns.CRS))

# add halite boundary - use buffering to smooth the outline
shape.buffer(1000).buffer(-1000).boundary.plot(
    ax=ax, color="black", linewidth=2
)

# wind farms
colours = ["lime", "lime", "black", "deepskyblue"]
for index, colour in zip(range(len(wind_farms)), colours):
    wind_farms.iloc[[index]].to_crs(fns.CRS).plot(
        ax=ax,
        hatch="///",
        facecolor="none",
        edgecolor=colour,
        linewidth=2,
        zorder=2,
    )
legend_handles = [
    mpatches.Patch(
        facecolor="none",
        hatch="////",
        edgecolor=colours[1:][x],
        label=list(wind_farms.dissolve("Name_")["Name"])[x],
    )
    for x in range(len(wind_farms.dissolve("Name_")))
]

legend_handles.append(
    mpatches.Patch(
        facecolor="none",
        edgecolor="black",
        label="Kish Basin halite",
        linewidth=2,
    )
)

weibull_k.to_crs(fns.CRS).plot(
    column="Value", cmap="flare", figsize=(6, 6), legend=True, ax=ax, zorder=1
)

cx.add_basemap(ax, crs=fns.CRS, source=cx.providers.CartoDB.Voyager, zoom=10)
ax.gridlines(
    draw_labels={"bottom": "x", "left": "y"}, alpha=0.25, color="darkslategrey"
)
ax.add_artist(
    ScaleBar(1, box_alpha=0, location="lower right", color="darkslategrey")
)
ax.legend(handles=legend_handles, loc="upper right")

plt.title(None)
plt.tight_layout()
plt.show()

# areas intersecting with wind farms
weibull_c = weibull_c.overlay(wind_farms, how="intersection")
weibull_k = weibull_k.overlay(wind_farms, how="intersection")

# average c and k over wind farms
weibull_c = wind_farms.dissolve(by="Name_").merge(
    weibull_c.dissolve(by="Name_", aggfunc={"Value": ["min", "max", "mean"]}),
    on="Name_",
)
weibull_k = wind_farms.dissolve(by="Name_").merge(
    weibull_k.dissolve(by="Name_", aggfunc={"Value": ["min", "max", "mean"]}),
    on="Name_",
)

weibull_c[["Name", ("Value", "min"), ("Value", "max"), ("Value", "mean")]]

weibull_k[["Name", ("Value", "min"), ("Value", "max"), ("Value", "mean")]]

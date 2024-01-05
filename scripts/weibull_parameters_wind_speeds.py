#!/usr/bin/env python
# coding: utf-8

# # Weibull Parameters of Wind Speeds 2001 to 2010 - 150m above ground level
#
# <https://data.gov.ie/dataset/weibull-parameters-wind-speeds-2001-to-2010-150m-above-ground-level>

import os
from zipfile import ZipFile

import cartopy.crs as ccrs
import contextily as cx
import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_scalebar.scalebar import ScaleBar

from src import data as rd
from src import functions as fns

# base data download directory
DATA_DIR = os.path.join("data", "weibull-parameters-wind-speeds")

URL = (
    "https://seaiopendata.blob.core.windows.net/wind/"
    "Weibull_150m_params_ITM.zip"
)

FILE_NAME = "Weibull_150m_params_ITM.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

weibull_c = rd.read_shapefile_from_zip(
    data_path=os.path.join(DATA_FILE), endswith="c_ITM.shp"
)

weibull_k = rd.read_shapefile_from_zip(
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

ds, extent = rd.read_dat_file(dat_path=os.path.join("data", "kish-basin"))

xmin, ymin, xmax, ymax = extent.total_bounds

ax = weibull_c.to_crs(3857).plot(
    column="Value",
    cmap="flare",
    figsize=(6, 6),
    legend=True,
    legend_kwds={"label": "c"},
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
    legend_kwds={"label": "k"},
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# wind farms in the area
wind_farms = fns.constraint_wind_farm(
    data_path=os.path.join(
        "data", "wind-farms", "wind-farms-foreshore-process.zip"
    ),
    dat_extent=extent,
)

# combine Codling wind farm polygons
wind_farms["Name_"] = wind_farms["Name"].str.split(expand=True)[0]

# shape of the halite
shape = rd.halite_shape(dat_xr=ds)

# land boundary
land = rd.read_shapefile_from_zip(
    data_path=os.path.join(
        "data", "boundaries", "osi-provinces-ungeneralised-2019.zip"
    )
)

land = land.dissolve().to_crs(rd.CRS)

# crop to wind farm and basin extent
extent_wf = gpd.GeoDataFrame(
    geometry=(
        gpd.GeoDataFrame(geometry=wind_farms.dissolve().envelope)
        .overlay(gpd.GeoDataFrame(geometry=extent), how="union")
        .dissolve()
        .envelope
    )
)
weibull_c = weibull_c.to_crs(rd.CRS).overlay(extent_wf, how="intersection")
weibull_k = weibull_k.to_crs(rd.CRS).overlay(extent_wf, how="intersection")

# crop land boundary from c and k
weibull_c = weibull_c.overlay(land, how="difference")
weibull_k = weibull_k.overlay(land, how="difference")


def plot_map(df, label):
    plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.epsg(rd.CRS))

    # add halite boundary - use buffering to smooth the outline
    shape.buffer(1000).buffer(-1000).boundary.plot(
        ax=ax, color="black", linewidth=2
    )

    # wind farms
    colours = ["lime", "lime", "darkslategrey", "deepskyblue"]
    for index, colour in zip(range(len(wind_farms)), colours):
        wind_farms.iloc[[index]].plot(
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

    df.plot(
        column="Value",
        cmap="flare",
        figsize=(6, 6),
        legend=True,
        ax=ax,
        zorder=1,
        legend_kwds={"label": label},
    )

    cx.add_basemap(
        ax, crs=rd.CRS, source=cx.providers.CartoDB.Voyager, zoom=10
    )
    ax.gridlines(
        draw_labels={"bottom": "x", "left": "y"},
        alpha=0.25,
        color="darkslategrey",
    )
    ax.add_artist(
        ScaleBar(1, box_alpha=0, location="lower right", color="darkslategrey")
    )
    ax.legend(handles=legend_handles, loc="upper right")

    plt.title(None)
    plt.tight_layout()
    plt.show()


plot_map(weibull_c, "c")

plot_map(weibull_k, "k")

# areas intersecting with wind farms
weibull_c = weibull_c.overlay(wind_farms, how="intersection")
weibull_k = weibull_k.overlay(wind_farms, how="intersection")

# compute c and k over wind farms
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

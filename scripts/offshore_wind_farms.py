#!/usr/bin/env python
# coding: utf-8

# # Offshore wind farms
#
# - <https://data.gov.ie/dataset/wind-farms-foreshore-process>
# - <https://data-housinggovie.opendata.arcgis.com/maps/housinggovie::wind-farms-foreshore-process>

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
    "803a4ecc22aa4cc09111072a0bbc4fac_2/downloads/"
    "data?format=shp&spatialRefId=4326&where=1%3D1"
)

FILE_NAME = "wind-farms-foreshore-process.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

plt.rcParams["xtick.major.size"] = 0
plt.rcParams["ytick.major.size"] = 0
plt.rcParams["xtick.minor.size"] = 0
plt.rcParams["ytick.minor.size"] = 0

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

wind_farms = rd.read_shapefile_from_zip(data_path=os.path.join(DATA_FILE))

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

ds, extent = rd.read_dat_file(dat_path=DATA_DIR)

# use extent bounds
xmin, ymin, xmax, ymax = extent.total_bounds

# shape of the halite
shape = rd.halite_shape(dat_xr=ds)

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

# combine Codling wind farm polygons
wind_farms_["Name_"] = wind_farms_["Name"].str.split(expand=True)[0]

wind_farms_ = wind_farms_.dissolve(by="Name_")

# remove abbreviation from name
wind_farms_.at["North", "Name"] = wind_farms_.at["North", "Name"].split(" (")[
    0
]

plt.figure(figsize=(11, 11))
ax = plt.axes(projection=ccrs.epsg(rd.CRS))

ds.max(dim="halite")["Thickness"].plot.contourf(
    cmap="jet",
    alpha=0.65,
    robust=True,
    levels=15,
    cbar_kwargs={
        "label": "Maximum Halite Thickness [m]",
        "aspect": 25,
        "pad": 0.035,
    },
)

plt.xlim(xmin - 8550, xmax + 1000)
# plt.ylim(ymin - 10500, ymax + 10500)

# wind farms
colours = ["firebrick", "black", "royalblue"]
legend_handles = []
for index, colour in zip(range(len(wind_farms_)), colours):
    wind_farms_.iloc[[index]].to_crs(rd.CRS).to_crs(rd.CRS).plot(
        ax=ax, hatch="///", facecolor="none", linewidth=2, edgecolor=colour
    )
    legend_handles.append(
        mpatches.Patch(
            facecolor="none",
            hatch="///",
            edgecolor=colour,
            label=wind_farms_.iloc[[index]]["Name"].values[0],
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
for index, colour in zip(range(len(wind_farms_)), colours):
    wind_farms_.iloc[[index]].to_crs(rd.CRS).to_crs(rd.CRS).plot(
        ax=ax, hatch="///", facecolor="none", linewidth=2, edgecolor=colour
    )
    legend_handles.append(
        mpatches.Patch(
            facecolor="none",
            hatch="///",
            edgecolor=colour,
            label=wind_farms_.iloc[[index]]["Name"].values[0],
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
wind_farms_ = wind_farms_.dissolve("Name_").reset_index().to_crs(rd.CRS)

for i in range(len(wind_farms_)):
    print(
        wind_farms_.iloc[[i]]["Name"].values[0],
        "is",
        wind_farms_.iloc[[i]]
        .distance(shape["geometry"], align=False)
        .values[0],
        "m away from Kish Bank",
    )

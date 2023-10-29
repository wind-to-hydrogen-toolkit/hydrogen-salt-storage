#!/usr/bin/env python
# coding: utf-8

# # Wind Farms (Foreshore Process)
#
# <https://data.gov.ie/dataset/wind-farms-foreshore-process>

import glob
import os
from datetime import datetime, timezone
from zipfile import ZipFile

import cartopy.crs as ccrs
import contextily as cx
import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import pooch
import shapely
from geocube.api.core import make_geocube
from matplotlib_scalebar.scalebar import ScaleBar

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

crs = 23029


def read_dat_file(dat_path: str, dat_crs):
    """
    Read XYZ data layers into an Xarray dataset
    """

    gdf = {}
    for dat_file in glob.glob(os.path.join(dat_path, "*.dat")):
        # read each layer as individual dataframes into a dictionary
        gdf[os.path.split(dat_file)[1][:-4]] = pd.read_fwf(
            dat_file, header=None, names=["X", "Y", "Z"]
        )

        # assign layer name to a column
        gdf[os.path.split(dat_file)[1][:-4]]["data"] = os.path.split(dat_file)[
            1
        ][:-4]

    # find data resolution
    gdf_xr = (
        gdf[os.path.split(dat_file)[1][:-4]].set_index(["X", "Y"]).to_xarray()
    )
    resx = gdf_xr["X"][1] - gdf_xr["X"][0]
    resy = gdf_xr["Y"][1] - gdf_xr["Y"][0]

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
    ds = make_geocube(
        vector_data=gdf,
        resolution=(-abs(resy), abs(resx)),
        align=(abs(resy / 2), abs(resx / 2)),
        group_by="data",
    )

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

    return ds, extent


ds, extent = read_dat_file(DATA_DIR, dat_crs=crs)

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
ax = plt.axes(projection=ccrs.epsg(crs))

# halite
ds.sel(data=[x for x in ds["data"].values if "Thickness XYZ" in x]).max(
    dim="data"
)["Z"].plot.contourf(
    cmap="jet",
    alpha=0.5,
    robust=True,
    levels=15,
    xlim=(xmin, xmax),
    ylim=(ymin - 4000, ymax + 4000),
    cbar_kwargs={"label": "Maximum Halite Thickness (m)"},
)

# wind farms
colours = ["firebrick", "forestgreen", "black", "royalblue"]
for index, colour in zip(range(len(wind_farms_)), colours):
    wind_farms_.iloc[[index]].to_crs(crs).plot(
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

cx.add_basemap(ax, crs=crs, source=cx.providers.CartoDB.Voyager)
ax.gridlines(
    draw_labels={"bottom": "x", "left": "y"}, alpha=0.25, color="darkslategrey"
)
ax.add_artist(
    ScaleBar(1, box_alpha=0, location="lower right", color="darkslategrey")
)
ax.legend(handles=legend_handles, loc="lower right", bbox_to_anchor=(1, 0.05))

plt.title("Wind Farms near Kish Basin")
# plt.title(None)
plt.tight_layout()
plt.show()

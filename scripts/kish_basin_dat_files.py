#!/usr/bin/env python
# coding: utf-8

# # Kish Basin Halite Data
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

# base data download directory
DATA_DIR = os.path.join("data", "kish-basin")
os.makedirs(DATA_DIR, exist_ok=True)

URL = "https://hyss.ie/wp-content/uploads/2023/07/Kish-Basin-dat-files.zip"
KNOWN_HASH = None
FILE_NAME = "Kish-Basin-dat-files.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# boundary data
ie = gpd.read_file(
    os.path.join("data", "boundaries.gpkg"), layer="NUTS_RG_01M_2021_4326_IE"
)

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

ZipFile(DATA_FILE).namelist()

# extract archive
try:
    z = ZipFile(DATA_FILE)
    z.extractall(DATA_DIR)
except BadZipFile:
    print("There were issues with the file", DATA_FILE)

# ## Map extent

with open(os.path.join(DATA_DIR, "Kish GIS Map Extent - Square.csv")) as f:
    print(f.read())

crs = 23029

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

extent.crs

ax = ie.to_crs(crs).plot(
    color="navajowhite",
    figsize=(7, 7),
    edgecolor="darkslategrey",
    linewidth=0.4,
)
extent.boundary.plot(ax=ax)

plt.title("Kish GIS Map Extent")
plt.text(
    500000,
    5.68e6,
    "© EuroGeographics for the administrative boundaries\n" "© HYSS",
)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# ## XYZ data


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

ds

ds.rio.crs


def plot_maps(plot_data):
    """
    Create facet contour plots of the Xarray dataset
    """

    fig = plot_data["Z"].plot.contourf(
        # x="x", y="y",
        col="data",
        cmap="jet",
        col_wrap=2,
        robust=True,
        levels=15,
        subplot_kws={"projection": ccrs.epsg(crs)},
        transform=ccrs.epsg(crs),
        xlim=(687000, 742000),
        ylim=(5888000, 5937000),
        cbar_kwargs={"aspect": 20, "pad": 0.02},
    )
    for axis in fig.axs.flat:
        ie.to_crs(crs).boundary.plot(
            ax=axis, edgecolor="darkslategrey", linewidth=0.5
        )
    # fig.set_titles("{value}", fontsize=10)
    # assign titles this way to prevent truncation/overflow
    for ax, title in zip(fig.axs.flat, plot_data["data"].values):
        # ax.set_title(title, fontsize=10, wrap=True)
        ax.set_title("\n".join(wrap(title, 34)), fontsize=10)
    plt.show()


# ### Halite thickness

plot_maps(ds.sel(data=[x for x in ds["data"].values if "Thickness XYZ" in x]))

# ### Halite thickness - zone of interest

plot_maps(ds.sel(data=[x for x in ds["data"].values if "Zone" in x]))

# ### Halite base depth

plot_maps(ds.sel(data=[x for x in ds["data"].values if "Base" in x]))

# ### Halite top depth

plot_maps(ds.sel(data=[x for x in ds["data"].values if "Top Depth" in x]))

# ### Halite top TWT (two-way thickness)

plot_maps(ds.sel(data=[x for x in ds["data"].values if "Millisecond" in x]))

# to create a 3D plot (triangular surface)
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# ds = ds.isel(data=0)
# ax.plot_trisurf(ds.X, ds.Y, ds.Z, cmap="Spectral_r")
# plt.axis("equal")
# plt.tight_layout()
# plt.show()

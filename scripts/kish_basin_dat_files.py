#!/usr/bin/env python
# coding: utf-8

# # Kish Bank Halite
#
# https://hyss.ie/

import os
from zipfile import BadZipFile, ZipFile
import geopandas as gpd
import matplotlib.pyplot as plt
import rioxarray as rxr
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import pooch
import glob
from datetime import datetime, timezone
from shapely.geometry import Polygon

# base data download directory
DATA_DIR = os.path.join("data", "kish-basin")
os.makedirs(DATA_DIR, exist_ok=True)

URL = "https://hyss.ie/wp-content/uploads/2023/07/Kish-Basin-dat-files.zip"
KNOWN_HASH = None
FILE_NAME = "Kish-Basin-dat-files.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# boundary data
ie = gpd.read_file(
    os.path.join("data", "boundaries", "ref-nuts-2021-01m.gpkg"),
    layer="NUTS_RG_01M_2021_4326_IE",
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

try:
    z = ZipFile(DATA_FILE)
    z.extractall(DATA_DIR)
except BadZipFile:
    print("There were issues with the file", DATA_FILE)

# ## Map extent

with open(os.path.join(DATA_DIR, "Kish GIS Map Extent - Square.csv")) as f:
    print(f.read())

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
    crs=23029,
)

# extent = pd.read_csv(
#     os.path.join(DATA_DIR, "Kish GIS Map Extent - Square.csv"), skiprows=2
# )
# extent["wkt"] = (
#     "POINT (" + extent[" X"].astype(str) + " " +
#     extent[" Y"].astype(str) + ")"
# )
# extent = gpd.GeoDataFrame(
#     extent, geometry=gpd.GeoSeries.from_wkt(extent["wkt"]), crs=23029
# )
# extent.drop(columns=["wkt"], inplace=True)

ax = ie.to_crs(23029).plot(
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


def read_dat_file(dat_path: str, crs: int = 23029):
    ds = {}
    for dat_file in glob.glob(dat_path):
        dat = pd.read_fwf(dat_file, header=None, names=["X", "Y", "Z"])
        dat["wkt"] = (
            "POINT ("
            + dat["X"].astype(str)
            + " "
            + dat["Y"].astype(str)
            + " "
            + dat["Z"].astype(str)
            + ")"
        )
        dat = gpd.GeoDataFrame(
            dat, geometry=gpd.GeoSeries.from_wkt(dat["wkt"]), crs=crs
        )
        dat.drop(columns=["wkt"], inplace=True)

        ds[os.path.split(dat_file)[1][:-4]] = (
            dat[["X", "Y", "Z"]].set_index(["X", "Y"]).to_xarray()
        )
        ds[os.path.split(dat_file)[1][:-4]] = (
            ds[os.path.split(dat_file)[1][:-4]]
            .assign_coords(
                data=(
                    os.path.split(dat_file)[1][:-4]
                    .replace(" XYZ ", "\n")
                    .replace(" -\n", "\n")
                )
            )
            .expand_dims(dim="data")
        )

    ds = xr.combine_by_coords(ds.values(), combine_attrs="override")
    ds.rio.write_crs(crs, inplace=True)

    return ds


ds = read_dat_file(os.path.join(DATA_DIR, "*.dat"))

ds


def plot_maps(plot_data):
    fig = plot_data["Z"].plot.contourf(
        x="X",
        y="Y",
        col="data",
        cmap="jet",
        col_wrap=2,
        robust=True,
        levels=15,
        subplot_kws={"projection": ccrs.epsg(23029)},
        transform=ccrs.epsg(23029),
        xlim=(687000, 742000),
        ylim=(5888000, 5937000),
        cbar_kwargs={"aspect": 20, "pad": 0.02},
    )
    for axis in fig.axs.flat:
        ie.to_crs(23029).boundary.plot(
            ax=axis, edgecolor="darkslategrey", linewidth=0.5
        )
    # fig.set_titles("{value}", fontsize=10)
    for ax, title in zip(fig.axs.flat, plot_data["data"].values):
        ax.set_title(title, fontsize=10)
    plt.show()


# ### Halite thickness

plot_maps(
    ds.sel(data=[x for x in ds["data"].values if "Thickness\nMeters" in x])
)

# ### Halite thickness - zone of interest

plot_maps(ds.sel(data=[x for x in ds["data"].values if "Zone" in x]))

# ### Halite base depth

plot_maps(ds.sel(data=[x for x in ds["data"].values if "Base" in x]))

# ### Halite top depth

plot_maps(ds.sel(data=[x for x in ds["data"].values if "Top Depth" in x]))

# ### Halite top TWT (two-way thickness)

plot_maps(ds.sel(data=[x for x in ds["data"].values if "Millisecond" in x]))

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# ds = ds.isel(data=0)
# ax.plot_trisurf(ds.X, ds.Y, ds.Z, cmap="Spectral_r")
# plt.axis("equal")
# plt.tight_layout()
# plt.show()

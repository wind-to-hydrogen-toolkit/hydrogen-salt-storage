#!/usr/bin/env python
# coding: utf-8

# # KIS-ORCA subsea cables
#
# <https://kis-orca.org/>
#
# **IMPORTANT:** There may be some incorrect name assignments as a simple
# backfill and multiline conversion was used and no further checks were done

import os
from datetime import datetime, timezone
from zipfile import BadZipFile, ZipFile

import contextily as cx
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pooch

# base data download directory
DATA_DIR = os.path.join("data", "kis-orca")
os.makedirs(DATA_DIR, exist_ok=True)

URL = "https://kis-orca.org/wp-content/uploads/2020/12/Olex_KIS-ORCA-v2023.zip"
KNOWN_HASH = None
FILE_NAME = "Olex_KIS-ORCA-v2023.zip"

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

# ## Read data

# extract the archive
try:
    z = ZipFile(DATA_FILE)
    z.extractall(DATA_DIR)
except BadZipFile:
    print("There were issues with the file", DATA_FILE)

DATA_FILE = os.path.join(DATA_DIR, ZipFile(DATA_FILE).namelist()[0])

# extract gz file using gzip
# https://www.gnu.org/software/gzip/
os.system(f"gzip -d < {DATA_FILE} > {DATA_FILE[:-3]}")

DATA_FILE = DATA_FILE[:-3]

for n, line in enumerate(open(DATA_FILE, "r", encoding="ISO-8859-15")):
    if n < 11:
        print(line[:-1])

# coordinates
with open(f"{DATA_FILE}_data.txt", "w", encoding="utf-8") as outfile:
    for n, line in enumerate(open(DATA_FILE, "r", encoding="ISO-8859-15")):
        if line[0].isdigit():
            outfile.write(f"{n} {line}")

# names
with open(f"{DATA_FILE}_names.txt", "w", encoding="utf-8") as outfile:
    for n, line in enumerate(open(DATA_FILE, "r", encoding="ISO-8859-15")):
        if "MTekst 2" in str(line):
            outfile.write(f"{n}, {line[18:]}")

data = pd.read_csv(
    f"{DATA_FILE}_data.txt",
    header=None,
    sep=" ",
    names=["y", "x", "z", "text"],
)

names = pd.read_csv(
    f"{DATA_FILE}_names.txt", header=None, names=["index", "name"]
)

data.head()

data.shape

names.head()

names.shape

names.set_index("index", inplace=True)

# ## Merge names with coordinates

names = names.reindex(range(max(data.index) + 1))

# handle missing data with backfill / forward fill
names = names.bfill()
names = names.ffill()

names.head()

names.shape

# merge names and data
data = pd.merge(data, names, left_index=True, right_index=True)

data.head()

data.shape

# ## Convert to geodataframe

# drop duplicate entries using coordinates
data = data.drop_duplicates(["y", "x"])

data.shape

# convert coords from minutes to degrees
# https://gis.stackexchange.com/a/241922
data["x"] = data["x"] / 60
data["y"] = data["y"] / 60

data.head()

# convert to geodataframe
data = gpd.GeoDataFrame(
    data, geometry=gpd.points_from_xy(data.x, data.y, crs=4326)
)

data.head()

data.crs

ax = data.to_crs(3857).plot(marker=".", markersize=1, figsize=(9, 9))
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# ## Dissolve geometries

# dissolve by name
data = data.dissolve("name").reset_index()

data.head()

data.shape

# ## Data for the Irish Sea

# extent
mask = (-7, 53, -4.5, 54)

data_ie = data.clip(mask)

data_ie

data_ie.shape

# ## Convert multipoint geometries to multilines

# convert all multi points to multi lines
data_ie_1 = data_ie[
    data_ie["geometry"].astype(str).str.contains("MULTI")
].reset_index(drop=True)

data_ie_1["geometry"] = gpd.GeoSeries.from_wkt(
    "LINESTRING ("
    + data_ie_1["geometry"].astype(str).str.split("(", expand=True)[1],
    crs=4326,
)

# merge multilines with remaining geometry data
data_ie = pd.concat(
    [
        data_ie[~data_ie["geometry"].astype(str).str.contains("MULTI")],
        data_ie_1,
    ]
).reset_index(drop=True)

data_ie

# get bounds of Irish Sea data
xmin, ymin, xmax, ymax = data_ie.to_crs(23029).total_bounds

# view plotter points in the Irish Sea
ax = data.to_crs(23029).plot(alpha=0.5, figsize=(9, 9))
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, crs=23029)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# remove obvious incorrect data lines - Hibernia Atlantic
data_ie = data_ie.drop([6]).reset_index(drop=True)

ax = (
    gpd.GeoDataFrame(geometry=data_ie.to_crs(23029).buffer(750))
    .dissolve()
    .plot(figsize=(12, 12), alpha=0.25, color="slategrey")
)
data_ie.to_crs(23029).plot(
    column="name",
    legend=True,
    ax=ax,
    cmap="jet",
    legend_kwds={"loc": "upper right"},
)
data.to_crs(23029).plot(
    ax=ax, marker="x", color="black", markersize=20, alpha=0.5
)
plt.xlim(xmin - 10000, xmax + 50000)
plt.ylim(ymin, ymax)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, crs=23029)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

data_ie.to_file(os.path.join(DATA_DIR, "KIS-ORCA.gpkg"))

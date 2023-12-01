#!/usr/bin/env python
# coding: utf-8

# # OREDP Boundaries

import os
from datetime import datetime, timezone
from zipfile import ZipFile

import contextily as cx
import geopandas as gpd
import matplotlib.pyplot as plt
import pooch
import seaborn as sns

# ## OREDP Assessment Zone
# 
# <https://www.isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.2212>

# base data download directory
DATA_DIR = os.path.join("data", "oredp-zone")
os.makedirs(DATA_DIR, exist_ok=True)

URL = (
    "https://atlas.marine.ie/midata/EnergyResourcesTidal/"
    "OREDP_Assessment_Zone.shapezip.zip"
)
KNOWN_HASH = None
FILE_NAME = "OREDP_Assessment_Zone.shapezip.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

# download data if necessary
if not os.path.isfile(DATA_FILE):
    pooch.retrieve(
        url=URL, known_hash=KNOWN_HASH, fname=FILE_NAME, path=DATA_DIR
    )

    with open(f"{DATA_FILE[:-13]}.txt", "w", encoding="utf-8") as outfile:
        outfile.write(
            f"Data downloaded on: {datetime.now(tz=timezone.utc)}\n"
            f"Download URL: {URL}"
        )

with open(f"{DATA_FILE[:-13]}.txt", encoding="utf-8") as f:
    print(f.read())

ZipFile(DATA_FILE).namelist()

data = gpd.read_file(
    os.path.join(
        f"zip://{DATA_FILE}!"
        + [x for x in ZipFile(DATA_FILE).namelist() if x.endswith(".shp")][0]
    )
)

data

data.shape

data.columns

data.crs

ax = data.to_crs(3857).plot(
    figsize=(7.5, 7.5),
    column="styleLayer",
    alpha=0.85,
    cmap=sns.color_palette("crest", as_cmap=True),
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
plt.title("OREDP Assessment Zone")

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# ## OREDP Study Area
# 
# <https://www.isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.2214>

URL = (
    "https://atlas.marine.ie/midata/EnergyResourcesTidal/"
    "OREDP_Study_Area.shapezip.zip"
)
FILE_NAME = "OREDP_Study_Area.shapezip.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# download data if necessary
if not os.path.isfile(DATA_FILE):
    pooch.retrieve(
        url=URL, known_hash=KNOWN_HASH, fname=FILE_NAME, path=DATA_DIR
    )

    with open(f"{DATA_FILE[:-13]}.txt", "w", encoding="utf-8") as outfile:
        outfile.write(
            f"Data downloaded on: {datetime.now(tz=timezone.utc)}\n"
            f"Download URL: {URL}"
        )

with open(f"{DATA_FILE[:-13]}.txt", encoding="utf-8") as f:
    print(f.read())

ZipFile(DATA_FILE).namelist()

data2 = gpd.read_file(
    os.path.join(
        f"zip://{DATA_FILE}!"
        + [x for x in ZipFile(DATA_FILE).namelist() if x.endswith(".shp")][0]
    )
)

data2

data2.shape

data2.columns

data2.crs

ax = data2.to_crs(3857).plot(
    figsize=(9, 9),
    column="localId",
    legend=True,
    legend_kwds={"loc": "lower right"},
    cmap="Accent",
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
plt.title("OREDP Study Area")

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()


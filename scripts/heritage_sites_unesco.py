#!/usr/bin/env python
# coding: utf-8

# # Heritage sites

import os
from datetime import datetime, timezone
from zipfile import ZipFile

import contextily as cx
import geopandas as gpd
import matplotlib.pyplot as plt
import pooch

# ## UNESCO Sites in Ireland
#
# <https://www.isde.ie/geonetwork/srv/eng/catalog.search#/metadata/69df8904-53df-4e1e-bddf-ab725a4060d4>

# base data download directory
DATA_DIR = os.path.join("data", "heritage")
os.makedirs(DATA_DIR, exist_ok=True)

URL = (
    "https://www.heritagecouncil.ie/content/files/"
    "UNESCO-Sites-in-Ireland-Shapefiles.zip"
)
KNOWN_HASH = None
FILE_NAME = "UNESCO-Sites-in-Ireland-Shapefiles.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

# download data if necessary
if not os.path.isfile(DATA_FILE):
    pooch.retrieve(
        url=URL, known_hash=KNOWN_HASH, fname=FILE_NAME, path=DATA_DIR
    )

    with open(f"{DATA_FILE[:-15]}.txt", "w", encoding="utf-8") as outfile:
        outfile.write(
            f"Data downloaded on: {datetime.now(tz=timezone.utc)}\n"
            f"Download URL: {URL}"
        )

with open(f"{DATA_FILE[:-15]}.txt") as f:
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

data.crs

ax = data.to_crs(3857).plot(
    column="Title",
    legend=True,
    cmap="tab20b",
    figsize=(9, 9),
    legend_kwds={"loc": "upper right"},
    linewidth=0.5,
    edgecolor="darkslategrey",
)
plt.xlim(-1.25e6, -0.2e6)
plt.ylim(6.6e6, 7.6e6)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=7)

plt.title("UNESCO Sites in Ireland")

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

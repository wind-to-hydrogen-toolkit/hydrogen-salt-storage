#!/usr/bin/env python
# coding: utf-8

# # Shipping

import os
from datetime import datetime, timezone
from zipfile import ZipFile

import contextily as cx
import geopandas as gpd
import matplotlib.pyplot as plt
import pooch

# ## Frequently Used Routes (300 gross tonnes and above)
#
# <https://data.gov.ie/dataset/frequently-used-routes-300-gross-tonnes-and-above1>

# base data download directory
DATA_DIR = os.path.join("data", "shipping")
os.makedirs(DATA_DIR, exist_ok=True)

URL = (
    "https://data-housinggovie.opendata.arcgis.com/datasets/"
    "housinggovie::frequently-used-routes-300-gross-tonnes-and-above.zip?"
    "outSR=%7B%22latestWkid%22%3A3857%2C%22wkid%22%3A102100%7D"
)
KNOWN_HASH = None
FILE_NAME = "shipping_frequently_used_routes.zip"

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

data = gpd.read_file(
    os.path.join(
        f"zip://{DATA_FILE}!"
        + [x for x in ZipFile(DATA_FILE).namelist() if x.endswith(".shp")][0]
    )
)

data.shape

data.head()

data.crs

data = data.sort_values("weight", ascending=False)

data.weight.unique()

data.key_.unique()

data["weight"] = data["weight"].astype(str)

data["legend"] = data.weight + ": " + data.key_

ax = data.to_crs(3857).plot(
    figsize=(7.5, 7.5),
    column="legend",
    legend=True,
    cmap="tab20",
    legend_kwds={"loc": "upper left"},
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=7)

plt.title("Frequently Used Routes (300 gross tonnes and above)")

plt.tick_params(labelbottom=False, labelleft=False)
# ax.legend(list(data.key_.unique()))
plt.tight_layout()
plt.show()

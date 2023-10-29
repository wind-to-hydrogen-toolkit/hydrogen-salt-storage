#!/usr/bin/env python
# coding: utf-8

# # Marine protected sites

import os
from datetime import datetime, timezone
from zipfile import ZipFile

import contextily as cx
import geopandas as gpd
import matplotlib.pyplot as plt
import pooch

# ## Dublin Bay Biosphere Marine Zones
#
# <https://data.gov.ie/dataset/dublin-bay-biosphere-marine-zones>

# base data download directory
DATA_DIR = os.path.join("data", "marine")
os.makedirs(DATA_DIR, exist_ok=True)

URL = (
    "https://data-housinggovie.opendata.arcgis.com/datasets/housinggovie::"
    "dublin-bay-biosphere-marine-zones.zip?"
    "outSR=%7B%22latestWkid%22%3A3857%2C%22wkid%22%3A102100%7D"
)
KNOWN_HASH = None
FILE_NAME = "dublin-bay-biosphere-marine-zones.zip"
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

data

data.shape

data.crs

ax = data.plot(
    figsize=(7.5, 7.5),
    column="zone",
    legend=True,
    cmap="Accent",
    alpha=0.55,
    legend_kwds={"loc": "upper left"},
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=11)
plt.title("Dublin Bay Biosphere")

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

ax = data.dissolve().plot(figsize=(7.5, 7.5), alpha=0.55)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=11)
plt.title("Dublin Bay Biosphere")

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

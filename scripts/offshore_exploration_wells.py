#!/usr/bin/env python
# coding: utf-8

# # Exploration Wells in the Irish Offshore
# 
# <https://www.isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.2171>

import os
from datetime import datetime, timezone
from zipfile import ZipFile

import contextily as cx
import geopandas as gpd
import matplotlib.pyplot as plt
import pooch

# base data download directory
DATA_DIR = os.path.join("data", "exploration-wells-irish-offshore")
os.makedirs(DATA_DIR, exist_ok=True)

URL = (
    "https://atlas.marine.ie/midata/EnergyResourcesExploration/"
    "Exploration_Wells_Irish_Offshore.shapezip.zip"
)
KNOWN_HASH = None
FILE_NAME = "Exploration_Wells_Irish_Offshore.shapezip.zip"

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

wells = gpd.read_file(
    os.path.join(
        f"zip://{DATA_FILE}!"
        + [x for x in ZipFile(DATA_FILE).namelist() if x.endswith(".shp")][0]
    )
)

wells.shape

wells.head()

wells.crs

ax = wells.to_crs(3857).plot(
    column="AREA",
    legend=True,
    cmap="tab20b",
    figsize=(7.5, 7.5),
    legend_kwds={"loc": "upper right"},
    linewidth=0.5,
    edgecolor="darkslategrey",
)
plt.xlim(-1.6e6, -0.2e6)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=6)

plt.title("Exploration Wells in the Irish Offshore (1970-2019)")

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()


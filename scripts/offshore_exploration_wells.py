#!/usr/bin/env python
# coding: utf-8

# # Exploration Wells in the Irish Offshore
#
# <https://www.isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.2171>

import os
from datetime import datetime, timezone
from zipfile import ZipFile

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

# boundary data
ie = gpd.read_file(
    os.path.join("data", "boundaries.gpkg"), layer="NUTS_RG_01M_2021_4326_IE"
)

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

ZipFile(DATA_FILE).namelist()

wells = gpd.read_file(
    os.path.join(f"zip://{DATA_FILE}!{FILE_NAME[:-4]}/{FILE_NAME[:-13]}.shp")
)

wells.shape

wells.head()

wells.crs

ax = ie.plot(
    color="navajowhite",
    figsize=(7.5, 7.5),
    edgecolor="darkslategrey",
    linewidth=0.4,
)
wells.plot(
    ax=ax,
    column="AREA",
    legend=True,
    cmap="tab20",
    marker=".",
    legend_kwds={"loc": "upper right", "bbox_to_anchor": (1.375, 0.725)},
)

plt.title("Exploration Wells in the Irish Offshore (1970-2019)")

plt.text(
    -10.5,
    49.4,
    "© EuroGeographics for the administrative boundaries\n"
    "© Petroleum Affairs Division",
)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

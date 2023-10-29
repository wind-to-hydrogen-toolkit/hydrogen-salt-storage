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

# ## UNESCO Global Geoparks and Biospheres
#
# <https://data.gov.ie/dataset/unesco-global-geoparks-and-biospheres>

# base data download directory
DATA_DIR = os.path.join("data", "heritage")
os.makedirs(DATA_DIR, exist_ok=True)

URL = (
    "https://data-housinggovie.opendata.arcgis.com/datasets/"
    "housinggovie::unesco-global-geoparks-and-biospheres.zip?"
    "outSR=%7B%22latestWkid%22%3A3857%2C%22wkid%22%3A102100%7D"
)
KNOWN_HASH = None
FILE_NAME = "unesco-global-geoparks-and-biospheres.zip"
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

data2 = gpd.read_file(
    os.path.join(
        f"zip://{DATA_FILE}!"
        + [x for x in ZipFile(DATA_FILE).namelist() if x.endswith(".shp")][0]
    )
)

data2

data2.shape

data2.crs

ax = data2.to_crs(3857).plot(
    figsize=(7.5, 7.5), column="Name", legend=True, cmap="Accent", alpha=0.85
)
plt.xlim(-1.2e6, -0.6e6)
plt.ylim(6.675e6, 7.55e6)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=7)
plt.title("UNESCO Global Geoparks and Biospheres")

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

# ## Shipwrecks in Irish Waters
#
# <https://isde.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.5131>

URL = (
    "https://gsi.geodata.gov.ie/downloads/Marine/Data/Downloads/Shapefiles/"
    "IE_GSI_MI_Shipwrecks_IE_Waters_WGS84_LAT.zip"
)
FILE_NAME = "IE_GSI_MI_Shipwrecks_IE_Waters_WGS84_LAT.zip"
DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

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

data3 = gpd.read_file(
    os.path.join(
        f"zip://{DATA_FILE}!"
        + [x for x in ZipFile(DATA_FILE).namelist() if x.endswith(".shp")][0]
    )
)

data3.head()

data3.shape

data3.crs

ax = data3.to_crs(3857).plot(
    figsize=(7.5, 7.5),
    edgecolor="darkslategrey",
    markersize=15,
    color="chartreuse",
)
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron)
plt.title("Shipwrecks in Irish Waters")

plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

#!/usr/bin/env python
# coding: utf-8

# # NUTS (Nomenclature of territorial units for statistics)
#
# <https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units/nuts>

import os
from datetime import datetime, timezone
from zipfile import BadZipFile, ZipFile

import geopandas as gpd
import matplotlib.pyplot as plt
import pooch

# base data download directory
DATA_DIR = os.path.join("data", "boundaries", "NUTS2021")
os.makedirs(DATA_DIR, exist_ok=True)

URL = (
    "https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/"
    "ref-nuts-2021-01m.shp.zip"
)
KNOWN_HASH = None
FILE_NAME = "ref-nuts-2021-01m.shp.zip"

# file name for the GeoPackage where the boundary vector layer will be saved
GPKG_BOUNDARY = os.path.join("data", "boundaries.gpkg")

DATA_DIR_TEMP = os.path.join(DATA_DIR, "temp")

os.makedirs(DATA_DIR_TEMP, exist_ok=True)

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# download data if necessary
if not os.path.isfile(DATA_FILE):
    pooch.retrieve(
        url=URL, known_hash=KNOWN_HASH, fname=FILE_NAME, path=DATA_DIR
    )

    with open(f"{DATA_FILE[:-8]}.txt", "w", encoding="utf-8") as outfile:
        outfile.write(
            f"Data downloaded on: {datetime.now(tz=timezone.utc)}\n"
            f"Download URL: {URL}"
        )

with open(f"{DATA_FILE[:-8]}.txt") as f:
    print(f.read())

ZipFile(DATA_FILE).namelist()

# extract the archive
try:
    z = ZipFile(DATA_FILE)
    z.extractall(DATA_DIR_TEMP)
except BadZipFile:
    print("There were issues with the file", DATA_FILE)

DATA_FILE = os.path.join(DATA_DIR_TEMP, "NUTS_RG_01M_2021_4326_LEVL_1.shp.zip")

ZipFile(DATA_FILE).namelist()

nuts = gpd.read_file(f"zip://{DATA_FILE}!NUTS_RG_01M_2021_4326_LEVL_1.shp")

nuts.head()

nuts.crs

nuts = nuts[nuts["NUTS_ID"].str.contains("IE0|UKN")]

nuts

nuts.total_bounds.round(2)

# ## Island of Ireland boundary

ie = nuts.copy()

ie = ie.dissolve(by="LEVL_CODE", as_index=False)

ie = ie[["geometry"]]

ie = ie.assign(NAME="Ireland")

DESCRIPTION = (
    "Boundary for the Island of Ireland generated using NUTS 2021 Level 1 "
    "boundaries"
)

ie = ie.assign(DESCRIPTION=DESCRIPTION)

ie

ie.plot(
    color="navajowhite",
    figsize=(7.5, 7.5),
    edgecolor="darkslategrey",
    linewidth=0.4,
)

plt.title("Boundary of the Island of Ireland")
plt.text(-8.75, 51.275, "© EuroGeographics for the administrative boundaries")
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

ie.to_file(GPKG_BOUNDARY, layer="NUTS_RG_01M_2021_4326_IE")

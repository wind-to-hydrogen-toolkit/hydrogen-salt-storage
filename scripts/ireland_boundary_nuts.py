#!/usr/bin/env python
# coding: utf-8

# # NUTS (Nomenclature of territorial units for statistics)
#
# <https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units/nuts>

# import libraries
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

# file name for the GeoPackage where the boundary vector layers will be saved
GPKG_BOUNDARY = os.path.join("data", "boundaries", "ref-nuts-2021-01m.gpkg")

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

ZipFile(DATA_FILE).namelist()

# extract the archive
try:
    z = ZipFile(DATA_FILE)
    z.extractall(DATA_DIR_TEMP)
except BadZipFile:
    print("There were issues with the file", DATA_FILE)

DATA_FILE = os.path.join(DATA_DIR_TEMP, "NUTS_RG_01M_2021_4326_LEVL_2.shp.zip")

nuts2 = gpd.read_file(f"zip://{DATA_FILE}!NUTS_RG_01M_2021_4326_LEVL_2.shp")

nuts2.head()

nuts2 = nuts2[nuts2["NUTS_ID"].str.contains("IE|UKN")]

nuts2

nuts2.total_bounds.round(2)

# ## Island of Ireland boundary

ie = nuts2.copy()

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

# ## Island of Ireland boundary in Irish transverse mercator
#
# Useful for plotting
#
# EPSG:2157
#
# See <https://www.gov.uk/government/publications/uk-geospatial-data-standards-register/national-geospatial-data-standards-register#standards-for-coordinate-reference-systems>

ie.to_crs(2157, inplace=True)

ie

ie.to_file(GPKG_BOUNDARY, layer="NUTS_RG_01M_2021_2157_IE")

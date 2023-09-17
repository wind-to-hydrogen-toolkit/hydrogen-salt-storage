#!/usr/bin/env python
# coding: utf-8

# # Wind Farms (Foreshore Process)
#
# https://data.gov.ie/dataset/wind-farms-foreshore-process

# import libraries
import os
from datetime import datetime, timezone
from zipfile import ZipFile
import geopandas as gpd
import matplotlib.pyplot as plt
import pooch

# base data download directory
DATA_DIR = os.path.join("data", "wind-farms-foreshore-process")
os.makedirs(DATA_DIR, exist_ok=True)

URL = (
    "https://opendata.arcgis.com/api/v3/datasets/"
    "803a4ecc22aa4cc09111072a0bbc4fac_2/downloads/"
    "data?format=shp&spatialRefId=4326&where=1%3D1"
)
KNOWN_HASH = None
FILE_NAME = "wind-farms-foreshore-process.zip"

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# boundary data
ie = gpd.read_file(
    os.path.join("data", "boundaries", "ref-nuts-2021-01m.gpkg"),
    layer="NUTS_RG_01M_2021_4326_IE",
)

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

ZipFile(DATA_FILE).namelist()

wind_farms = gpd.read_file(
    os.path.join(f"zip://{DATA_FILE}!Energy_Offshore_Renewable.shp")
)

wind_farms.crs

wind_farms.shape

wind_farms.columns

wind_farms[["Name", "Type", "MDM_Catego"]]

ax = ie.plot(
    color="navajowhite",
    figsize=(7.5, 7.5),
    edgecolor="darkslategrey",
    linewidth=0.4,
)
wind_farms.boundary.plot(ax=ax)

plt.title("Wind Farms (Foreshore Process)")
plt.text(
    -8.75,
    51.275,
    "© EuroGeographics for the administrative boundaries\n"
    "© Dept. of Housing, Local Government, and Heritage",
)
plt.tick_params(labelbottom=False, labelleft=False)
plt.tight_layout()
plt.show()

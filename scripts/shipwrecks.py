#!/usr/bin/env python
# coding: utf-8

# # Shipwrecks
#
# <https://data.marine.ie/geonetwork/srv/eng/catalog.search#/metadata/ie.marine.data:dataset.5131>

import os
from zipfile import ZipFile

import contextily as cx
import matplotlib.pyplot as plt

from h2ss import data as rd

plt.rcParams["xtick.major.size"] = 0
plt.rcParams["ytick.major.size"] = 0
plt.rcParams["xtick.minor.size"] = 0
plt.rcParams["ytick.minor.size"] = 0

# base data download directory
DATA_DIR = os.path.join("data", "shipwrecks")

FILE_NAME = "IE_GSI_MI_Shipwrecks_IE_Waters_WGS84_LAT.zip"

URL = (
    "https://gsi.geodata.gov.ie/downloads/Marine/Data/Downloads/Shapefiles/"
    + FILE_NAME
)

DATA_FILE = os.path.join(DATA_DIR, FILE_NAME)

# basemap cache directory
cx.set_cache_dir(os.path.join("data", "basemaps"))

rd.download_data(url=URL, data_dir=DATA_DIR, file_name=FILE_NAME)

ZipFile(DATA_FILE).namelist()

data = rd.read_shapefile_from_zip(data_path=os.path.join(DATA_FILE))

data.head()

data.shape

data.crs

ax = data.to_crs(3857).plot(
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
